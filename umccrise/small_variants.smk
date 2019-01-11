#################
#### Somatic ####
from os.path import isfile, join
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.reference_data import get_key_genes_set
from umccrise import package_path
import cyvcf2
import toml
import csv
from ngs_utils.file_utils import which
from vcf_stuff import iter_vcf


localrules: small_variants, prep_anno_toml, prep_giab_bed


# Call variants with 1%
# Apply to VarDict only: (INFO/QUAL * TUMOR_AF) >= 4
# Call ensemble
# First PoN round: remove PoN_CNT>2'
# Annotate with PCGR (VEP+known cancer databases)
#   Tier 1 - variants of strong clinical significance
#   Tier 2 - variants of potential clinical significance
#   Tier 3 - variants of unknown clinical significance
#   Tier 4 - other coding variants
#   Noncoding variants
# Tier 1-3 - keep all variants
# Tier 4 and noncoding - filter with:
#   Remove gnomad_AF <=0.02
#   Remove PoN_CNT>{0 if issnp else 1}'
#   Remove indels in "bad_promoter" tricky regions
#   Remove DP<25 & AF<5% in tricky regions: gc15, gc70to75, gc75to80, gc80to85, gc85, low_complexity_51to200bp, low_complexity_gt200bp, non-GIAB confident, unless coding in cancer genes


rule somatic_vcf_annotate:
    input:
        vcf = lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].name}-ensemble-annotated.vcf.gz')
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-ANNO.vcf.gz',
    params:
        genome = run.genome_build,
        work_dir = 'work/{batch}/small_variants/annotation'
    # group: "small_variants"
    shell:
        'anno_somatic_vcf {input.vcf} -g {params.genome} -o {output.vcf} -w {params.work_dir}'

rule somatic_vcf_filter:
    input:
        vcf = rules.somatic_vcf_annotate.output.vcf,
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-ANNO-FILT.vcf.gz',
    # group: "small_variants"
    shell:
        'filter_somatic_vcf {input.vcf} -o {output.vcf}'

rule somatic_vcf_filter_af10:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-ANNO-FILT-AF10.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic-ensemble-ANNO-FILT-AF10.vcf.gz.tbi',
    # group: "small_variants"
    shell:
        'bcftools filter -e "INFO/TUMOR_AF<0.1" --soft-filter af10 --mode + ' \
        '-Oz {input.vcf} -o {output.vcf} && tabix -f -p vcf -f {output.vcf}'

rule somatic_vcf_filter_pass:
    input:
        vcf = rules.somatic_vcf_filter_af10.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-ANNO-FILT-PASS.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic-ensemble-ANNO-FILT-PASS.vcf.gz.tbi',
    # group: "small_variants"
    shell:
        'bcftools view -f.,PASS {input.vcf} -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}'

##################
#### Germline ####
# Annotate any events found in ~200 cancer predisposition gene set.
rule germline_vcf_subset:  # {batch}
    input:
        vcf = lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].normal.name}{GERMLINE_SUFFIX}-ensemble-annotated.vcf.gz'),
    output:
        vcf = 'work/{batch}/small_variants/raw_normal-ensemble-predispose_genes.vcf.gz',
        tbi = 'work/{batch}/small_variants/raw_normal-ensemble-predispose_genes.vcf.gz.tbi',
    params:
        ungz = lambda wc, output: get_ungz_gz(output[0])[0]
    # group: "small_variants"
    run:
        pcgr_toml_fpath = join(package_path(), 'pcgr', 'cpsr.toml')
        genes = [g for g in toml.load(pcgr_toml_fpath)['cancer_predisposition_genes']]
        def func(rec):
            if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                return rec
        iter_vcf(input.vcf, output.vcf, func)

# Preparations: annotate TUMOR_X and NORMAL_X fields, remove non-standard chromosomes and mitochondria, remove non-PASSed calls.
# Suites for PCGR, but for all other processing steps too
rule germline_vcf_prep:
    input:
        vcf = rules.germline_vcf_subset.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-normal-ensemble-predispose_genes.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-normal-ensemble-predispose_genes.vcf.gz.tbi',
    # group: "small_variants"
    shell:
        'pcgr_prep {input.vcf} |'
        ' bcftools view -f.,PASS -Oz -o {output.vcf}'
        ' && tabix -f -p vcf {output.vcf}'


#############

rule small_variants:
    input:
        expand(rules.somatic_vcf_filter.output.vcf, batch=batch_by_name.keys()),
        expand(rules.somatic_vcf_filter_af10.output.vcf, batch=batch_by_name.keys()),
        expand(rules.somatic_vcf_filter_pass.output.vcf, batch=batch_by_name.keys()),
        expand(rules.germline_vcf_prep.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/small_variants.done'))