from os.path import join, abspath
from ngs_utils.logger import warn
from ngs_utils.reference_data import get_key_genes, get_key_genes_bed
from ngs_utils.file_utils import safe_mkdir
from reference_data import api as refdata

from umccrise import package_path


localrules: cancer_report, conda_list


# Subset to GiaB confident intervals
rule subset_to_giab:
    input:
        vcf = '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz',
    params:
        regions = refdata.get_ref_file(run.genome_build, key=['hmf_giab_conf'])
    output:
        'work/{batch}/cancer_report/afs/somatic-confident.vcf.gz'
    group: "rmd_prep"
    shell:
        'bcftools view {input.vcf} -R {params.regions} -Oz -o {output}'

# Split multiallelics to avoid R parsing issues
rule split_multiallelic:
    input:
        vcf = 'work/{batch}/cancer_report/afs/somatic-confident.vcf.gz',
        ref_fa = refdata.get_ref_file(run.genome_build, key='fa')
    output:
        'work/{batch}/cancer_report/afs/somatic-confident-singleallelic.vcf.gz'
    group: "rmd_prep"
    shell:
        'bcftools annotate -x ^INFO/TUMOR_AF {input.vcf} -Ob | '
        'bcftools norm -m \'-\' -f {input.ref_fa} -Ob | ' 
        'bcftools sort -Oz -o {output} '
        '&& tabix -p vcf {output}'

rule afs:
    input:
        'work/{batch}/cancer_report/afs/somatic-confident-singleallelic.vcf.gz'
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.rgid
    output:
        'work/{batch}/cancer_report/afs/af_tumor.txt'
    group: "rmd_prep"
    shell:
        'bcftools view {input} -s {params.tumor_name} -Ou | '
        '(printf "af\n"; bcftools query -f "%INFO/TUMOR_AF\\n") > {output} '
        '&& test -e {output}'

rule afs_keygenes:
    input:
        vcf = 'work/{batch}/cancer_report/afs/somatic-confident-singleallelic.vcf.gz',
        bed = get_key_genes_bed(run.genome_build, coding_only=True),
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.rgid
    output:
        'work/{batch}/cancer_report/afs/af_tumor_keygenes.txt'
    group: "rmd_prep"
    shell:
        'bcftools view -f .,PASS {input.vcf} -s {params.tumor_name} -Ov'
        ' | bedtools intersect -a stdin -b {input.bed} -header'
        ' | (printf "chrom\\tpos\\tid\\tref\\talt\\taf\\n" ; bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/TUMOR_AF\\n")'
        ' > {output} && test -e {output}'


rule conda_list:
    output:
        txt = 'work/conda_pkg_list.txt'
    params:
        env = ['""', '_pcgr', '_hmf', '_cancer_report'],
    shell:
        "for e in {params.env}; do conda list -p {env_path}$e "
        "| grep -v ^# "
        "| awk -v var=env$e '{{ print var, $0 }}' >> {output} ; done"


rule run_cancer_report:
    input:
        rmd_files_dir        = join(package_path(), 'rmd_files'),
        key_genes            = get_key_genes(),
        af_global            = rules.afs.output[0],
        af_keygenes          = rules.afs_keygenes.output[0],
        somatic_snv          = '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz',
        somatic_sv           = lambda wc: rules.prep_sv_tsv.output[0]
                               if (batch_by_name[wc.batch].sv_vcf and 'structural' in stages) else [],
        purple_cnv           = 'work/{batch}/purple/{batch}.purple.cnv.somatic.tsv',
        purple_gene_cnv      = 'work/{batch}/purple/{batch}.purple.cnv.gene.tsv',
        purple_purity        = 'work/{batch}/purple/{batch}.purple.purity.tsv',
        purple_qc            = 'work/{batch}/purple/{batch}.purple.qc',

        purple_circos_png    = 'work/{batch}/purple/plot/{batch}.circos.png',
        purple_input_png     = 'work/{batch}/purple/plot/{batch}.input.png',
        purple_purity_png    = 'work/{batch}/purple/plot/{batch}.purity.range.png',
        purple_segment_png   = 'work/{batch}/purple/plot/{batch}.segment.png',
        purple_clonality_png = 'work/{batch}/purple/plot/{batch}.somatic.clonality.png',
        purple_ploidy_png    = 'work/{batch}/purple/plot/{batch}.somatic.png',
        purple_rainfall_png  = 'work/{batch}/purple/plot/{batch}.somatic.rainfall.png',
        purple_baf_png       = 'work/{batch}/purple/circos_baf/{batch}.circos_baf.png',

        wait_for_integration_sites = get_integration_sites_tsv_fn,
        oncoviral_present_viruses = 'work/{batch}/oncoviruses/present_viruses.txt',

        conda_list           = rules.conda_list.output.txt,

    params:
        report_rmd = 'cancer_report.Rmd',
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.rgid,
        work_dir = os.getcwd(),
        output_file = lambda wc, output: join(os.getcwd(), output[0]),
        rmd_genome_build = 'hg19' if run.genome_build in ['GRCh37', 'hg19'] else run.genome_build,
        af_global       = lambda wc, input: abspath(input.af_global),
        af_keygenes     = lambda wc, input: abspath(input.af_keygenes),
        somatic_snv     = lambda wc, input: abspath(input.somatic_snv),
        somatic_sv      = lambda wc, input: abspath(input.somatic_sv) if input.somatic_sv else 'NA',
        purple_gene_cnv = lambda wc, input: abspath(input.purple_gene_cnv),
        purple_cnv      = lambda wc, input: abspath(input.purple_cnv),
        purple_purity   = lambda wc, input: abspath(input.purple_purity),
        purple_qc       = lambda wc, input: abspath(input.purple_qc),
        conda_list      = lambda wc, input: abspath(input.conda_list),
    output:
        report_html = '{batch}/{batch}_cancer_report.html',
        rmd_tmp_dir = directory('work/{batch}/cancer_report/rmd_files'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000
        # TODO: memory based on the mutation number. E.g. over 455k tumor mutations need over 5G
    run:
        oncoviral_present_viruses = abspath(input.oncoviral_present_viruses)
        oncoviral_breakpoints_tsv = ''
        with open(input.oncoviral_present_viruses) as f:
            if [v for v in f.read().strip().split(',') if v]:
                oncoviral_breakpoints_tsv = abspath(f'work/{wildcards.batch}/oncoviruses/oncoviral_breakpoints.tsv')

        shell('cp -r {input.rmd_files_dir} {output.rmd_tmp_dir}')
        shell('mkdir -p {output.rmd_tmp_dir}/img')
        for img_path in [
            input.purple_circos_png,
            input.purple_input_png,
            input.purple_purity_png,
            input.purple_segment_png,
            input.purple_clonality_png,
            input.purple_ploidy_png,
            input.purple_rainfall_png,
            input.purple_baf_png,
        ]:
            shell('cp ' + img_path + ' {output.rmd_tmp_dir}/img/')
        shell(conda_cmd.format('cancer_report') + """
cd {output.rmd_tmp_dir} && \
Rscript -e "rmarkdown::render('{params.report_rmd}', \
output_file='{params.output_file}', \
params=list( \
tumor_name='{params.tumor_name}', \
batch_name='{wildcards.batch}', \
genome_build='{params.rmd_genome_build}', \
key_genes='{input.key_genes}', \
af_global='{params.af_global}', \
af_keygenes='{params.af_keygenes}', \
somatic_snv='{params.somatic_snv}', \
somatic_sv='{params.somatic_sv}', \
purple_gene_cnv='{params.purple_gene_cnv}', \
purple_cnv='{params.purple_cnv}', \
purple_purity='{params.purple_purity}', \
purple_qc='{params.purple_qc}', \
oncoviral_present_viruses='{oncoviral_present_viruses}', \
oncoviral_breakpoints_tsv='{oncoviral_breakpoints_tsv}', \
conda_list='{params.conda_list}' \
))" ; \
cd {params.work_dir} ; \
""")


#############

rule cancer_report:
    input:
        expand(rules.run_cancer_report.output[0], batch=batch_by_name.keys())
    output:
        temp(touch('log/cancer_report.done'))

