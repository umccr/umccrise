#################
#### Somatic ####
from os.path import isfile, join
from ngs_utils.file_utils import get_ungz_gz
from umccrise import package_path
import cyvcf2
import toml


localrules: small_variants


# Preparations: annotate TUMOR_X and NORMAL_X fields, remove non-standard chromosomes and mitochondria, remove non-PASSed calls.
# Suites for PCGR, but for all other processing steps too
rule somatic_vcf_prep:  # {batch}
    input:
        vcf = lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].name}-ensemble-annotated.vcf.gz')
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-prep.vcf.gz',
        tbi = 'work/{batch}/small_variants/somatic-ensemble-prep.vcf.gz.tbi'
    shell:
        'pcgr_prep {input.vcf} |'
        ' bcftools view -f.,PASS -Oz -o {output.vcf}'
        ' && tabix -p vcf {output.vcf}'

# Bcbio doesn't properly filter Strelka2 and MuTect2 calls, so need to post-filter
rule somatic_vcf_filter_af:  # {batch}
    input:
        vcf = rules.somatic_vcf_prep.output.vcf
    params:
        min_af = 0.1
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-prep-min_af.vcf.gz',
        tbi = 'work/{batch}/small_variants/somatic-ensemble-prep-min_af.vcf.gz.tbi'
    shell:
        'bcftools filter -e "TUMOR_AF<{params.min_af}" {input.vcf} -Oz -o {output.vcf} && tabix {output.vcf}'

rule somatic_vcf_pon:  # {batch}
    input:
        vcf = rules.somatic_vcf_filter_af.output.vcf,
        tbi = rules.somatic_vcf_filter_af.output.tbi
    params:
        genome_build = run.genome_build,
        ht = 1
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-pon_softfiltered.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic-ensemble-pon_softfiltered.vcf.gz.tbi'
    run:
        if pon_dir:
            shell('pon_anno {input.vcf} -h {params.ht} -o {output.vcf} --pon-dir ' + pon_dir + ' '
                  '&& tabix -p vcf {output.vcf}')
        else:
            shell('cp {input.vcf} {output.vcf} && cp {input.tbi} {output.tbi}')

rule somatic_vcf_pon_pass:  # {batch}
    input:
        rules.somatic_vcf_pon.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-pon_hardfiltered.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic-ensemble-pon_hardfiltered.vcf.gz.tbi'
    shell:
        'bcftools view -f.,PASS -Oz {input} -o {output.vcf}'
        ' && tabix -p vcf {output.vcf}'

##################
#### Germline ####
# Annotate any events found in ~200 cancer predisposition gene set.
rule germline_vcf_subset:  # {batch}
    input:
        vcf = lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].normal.name}{GERMLINE_SUFFIX}-ensemble-annotated.vcf.gz'),
    output:
        vcf = 'work/{batch}/small_variants/raw_normal-ensemble-predispose_genes.vcf.gz',
        tbi = 'work/{batch}/small_variants/raw_normal-ensemble-predispose_genes.vcf.gz.tbi'
    params:
        ungz = lambda wc, output: get_ungz_gz(output[0])[0]
    run:
        pcgr_toml_fpath = join(package_path(), 'pcgr', 'cpsr.toml')
        genes = [g for g in toml.load(pcgr_toml_fpath)['cancer_predisposition_genes']]
        inp_vcf = cyvcf2.VCF(input.vcf)
        ungz = get_ungz_gz(output.vcf)[0]
        out_vcf = cyvcf2.Writer(ungz, inp_vcf)
        out_vcf.write_header()
        for rec in inp_vcf:
            if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                out_vcf.write_record(rec)
        out_vcf.close()
        shell('bgzip {ungz} && tabix -p vcf {output.vcf}')

# Preparations: annotate TUMOR_X and NORMAL_X fields, remove non-standard chromosomes and mitochondria, remove non-PASSed calls.
# Suites for PCGR, but for all other processing steps too
rule germline_vcf_prep:
    input:
        vcf = rules.germline_vcf_subset.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-normal-ensemble-predispose_genes.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-normal-ensemble-predispose_genes.vcf.gz.tbi'
    shell:
        'pcgr_prep {input.vcf} |'
        ' bcftools view -f.,PASS -Oz -o {output.vcf}'
        ' && tabix -p vcf {output.vcf}'


#############

rule small_variants:
    input:
        expand(rules.somatic_vcf_pon_pass.output, batch=batch_by_name.keys()),
        expand(rules.germline_vcf_prep.output, batch=batch_by_name.keys())
    output:
        temp(touch('log/small_variants.done'))