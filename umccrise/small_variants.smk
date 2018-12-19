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


# Preparations: annotate TUMOR_X and NORMAL_X fields for PCGR, remove non-standard chromosomes and mitochondria, remove non-PASSed calls
rule somatic_vcf_prep:
    input:
        vcf = lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].name}-ensemble-annotated.vcf.gz')
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-prep.vcf.gz',
        tbi = 'work/{batch}/small_variants/somatic-ensemble-prep.vcf.gz.tbi'
    # group: "small_variants"
    shell:
        'pcgr_prep {input.vcf}'
        ' | bcftools view -f.,PASS'
        ' -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'

rule somatic_vcf_pon_anno:
    input:
        vcf = rules.somatic_vcf_prep.output.vcf,
        tbi = rules.somatic_vcf_prep.output.tbi,
    params:
        genome_build = run.genome_build,
        pon_hits = 3,
    output:
        vcf = 'work/{batch}/small_variants/pon_anno/somatic-ensemble.vcf.gz',
    # group: "small_variants"
    shell:
        'pon_anno {input.vcf} --pon-dir ' + pon_dir + ' | bgzip -c > {output.vcf} && tabix -p vcf {output.vcf}'
        # ' | bcftools filter -e "INFO/PoN_CNT>={params.pon_hits}" --soft-filter PoN --mode + -Oz -o {output.vcf}' \
        # ' && tabix -p vcf {output.vcf} '

rule somatic_vcf_keygenes:
    input:
        vcf = rules.somatic_vcf_pon_anno.output.vcf,
    output:
        vcf = 'work/{batch}/small_variants/keygenes/somatic-ensemble.vcf.gz',
    # group: "small_variants"
    run:
        genes = get_key_genes_set()
        def func(rec):
            if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                return rec
        iter_vcf(input.vcf, output.vcf, func)

rule somatic_vcf_pcgr_ready:
    input:
        full_vcf = rules.somatic_vcf_pon_anno.output.vcf,
        keygenes_vcf = rules.somatic_vcf_keygenes.output.vcf,
    output:
        vcf = 'work/{batch}/small_variants/pcgr_input/somatic-ensemble.vcf.gz',
    # group: "small_variants"
    run:
        total_vars = int(subprocess.check_output(f'bcftools view -H {input.full_vcf} | wc -l', shell=True).strip())
        vcf = input.full_vcf if total_vars <= 500_000 else input.keygenes_vcf  # to avoid PCGR choking on too many variants
        shell(f'bcftools annotate -x INFO/ANN {vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}')

rule somatic_vcf_pcgr_round1:
    input:
        vcf = rules.somatic_vcf_pcgr_ready.output.vcf,
        pcgr_data = pcgr_data,
    output:
        tiers = 'work/{batch}/small_variants/pcgr_round1/{batch}-somatic.pcgr.snvs_indels.tiers.tsv',
    params:
        output_dir = 'work/{batch}/small_variants/pcgr_round1',
        genome = run.genome_build,
        sample_name = '{batch}-somatic',
        opt='--no-docker' if not which('docker') else '',
    # group: "small_variants"
    resources:
        mem_mb = 20000
    shell:
        conda_cmd.format('pcgr') + which('pcgr') +
        ' {input.vcf} -g {params.genome} -o {params.output_dir} -s {params.sample_name} '
        '{params.opt} --pcgr-data {input.pcgr_data}'

rule somatic_vcf_pcgr_anno:
    input:
        vcf = rules.somatic_vcf_pcgr_ready.output.vcf,
        tiers = rules.somatic_vcf_pcgr_round1.output.tiers,
    output:
        vcf = 'work/{batch}/small_variants/annotation/pcgr/somatic-ensemble.vcf.gz',
    # group: "small_variants"
    run:
        gene_by_snp = dict()
        tier_by_snp = dict()
        with open(input.tiers) as f:
            reader = csv.DictReader(f, delimiter='\t', fieldnames=f.readline().strip().split('\t'))
            for row in reader:
                k = row['GENOMIC_CHANGE']
                gene_by_snp[k] = row['SYMBOL']
                tier_by_snp[k] = row['TIER']
        def func_hdr(vcf):
            vcf.add_info_to_header({'ID': 'PCGR_TIER', 'Description': 'TIER as reported by PCGR in .snvs_indels.tiers.tsv file', 'Type': 'String', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_GENE', 'Description': 'Gene symbol as reported by PCGR in .snvs_indels.tiers.tsv file', 'Type': 'String', 'Number': '1'})
        def func(rec):
            k = f'{rec.CHROM}:g.{rec.POS}{rec.REF}>{rec.ALT[0]}'
            rec.INFO['PCGR_TIER'] = tier_by_snp.get(k, '.').replace(' ', '_')
            rec.INFO['PCGR_GENE'] = gene_by_snp.get(k, '.')
            return rec
        iter_vcf(input.vcf, output.vcf, func, func_hdr)

rule prep_giab_bed:
    input:
        get_ref_file(run.genome_build, ['truth_sets', 'giab', 'bed'])
    output:
        'work/giab_conf.bed.gz'
    shell:
        'cat {input} | bgzip -c > {output} && tabix -p bed {output}'

rule prep_hmf_hotspots:
    input:
        get_ref_file(run.genome_build, key='hmf_hotspot'),
    output:
        'work/hmf_hotspot.vcf.gz'
    params:
        ungz = 'work/hmf_hotspot.vcf'
    shell: """
echo "##fileformat=VCFv4.2" >> {params.ungz} && 
echo "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO" >> {params.ungz} && 
gunzip -c {input} | py -x "print('\t'.join([x.split()[0], x.split()[1], '.', x.split()[2], x.split()[3], '.', '.', 'HS']))" >> {params.ungz} &&
bgzip {params.ungz} && 
tabix -p vcf {output}
"""

rule prep_anno_toml:
    input:
        tricky_bed      = get_ref_file(run.genome_build, key='tricky'),
        giab_conf_bed   = rules.prep_giab_bed.output[0],
        gnomad_vcf      = get_ref_file(run.genome_build, key='gnomad'),
        hmf_hotspots    = rules.prep_hmf_hotspots.output[0],
        hmf_giab        = get_ref_file(run.genome_build, key='hmf_giab_conf'),
        hmf_mappability = get_ref_file(run.genome_build, key='hmf_mappability'),
    output:
        'work/tricky_vcfanno.toml'
    params:
        toml_text  = lambda wc, input, output: f'''
[[annotation]]
file = "{input.tricky_bed}"
names = ["TRICKY"]
columns = [4]
ops = ["self"]

[[annotation]]
file = "{input.giab_conf_bed}"
names = ["GIAB_CONF"]
columns = [3]
ops = ["flag"]

[[annotation]]
file="{input.gnomad_vcf}"
fields = ["AF"]
names = ["gnomAD_AF"]
ops = ["self"]

[[annotation]]
file = "{input.hmf_hotspots}"
fields = ["HS"]
names = ["HMF_HOTSPOT"]
ops = ["flag"]

[[annotation]]
file = "{input.hmf_giab}"
names = ["HMF_GIAB_CONF"]
columns = [3]
ops = ["flag"]

[[annotation]]
file = "{input.hmf_mappability}"
names = ["HMF_MAPPABILITY"]
columns = [5]
ops = ["self"]

'''.replace('\n', r'\\n').replace('"', r'\"'),
    shell:
        'printf "{params.toml_text}" > {output}'

rule somatic_vcf_regions_anno:
    input:
        vcf = rules.somatic_vcf_pcgr_anno.output.vcf,
        toml = rules.prep_anno_toml.output[0]
    output:
        vcf = 'work/{batch}/small_variants/annotation/regions/somatic-ensemble.vcf.gz',
    # group: "small_variants"
    shell:
        'vcfanno {input.toml} {input.vcf} | bgzip -c > {output} && tabix -p vcf -f {output}'

rule somatic_vcf_filter:
    input:
        vcf = rules.somatic_vcf_regions_anno.output.vcf,
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-filt.vcf.gz',
    # group: "small_variants"
    shell:
        'filter_somatic_vcf {input.vcf} -o {output.vcf}'

rule somatic_vcf_filter_pass:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-filt.pass.vcf.gz',
    # group: "small_variants"
    shell:
        'bcftools view -f.,PASS {input.vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'

##################
#### Germline ####
# Annotate any events found in ~200 cancer predisposition gene set.
rule germline_vcf_subset:  # {batch}
    input:
        vcf = lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].normal.name}{GERMLINE_SUFFIX}-ensemble-annotated.vcf.gz'),
    output:
        vcf = 'work/{batch}/small_variants/raw_normal-ensemble-predispose_genes.vcf.gz',
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
    # group: "small_variants"
    shell:
        'pcgr_prep {input.vcf} |'
        ' bcftools view -f.,PASS -Oz -o {output.vcf}'
        ' && tabix -p vcf {output.vcf}'


#############

rule small_variants:
    input:
        expand(rules.somatic_vcf_filter.output.vcf, batch=batch_by_name.keys()),
        expand(rules.somatic_vcf_filter_pass.output.vcf, batch=batch_by_name.keys()),
        expand(rules.germline_vcf_prep.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/small_variants.done'))