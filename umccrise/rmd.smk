from os.path import join

from ngs_utils.logger import warn
from ngs_utils.reference_data import get_key_genes, get_key_genes_bed
from ngs_utils.file_utils import safe_mkdir
import glob
from umccrise import package_path


localrules: rmd


## Allelic frequencies
# AF is not yet integrated into PCGR or CGI. We can extract those from VarDict or Mutect2 for plotting purposes,
# but still need to look up AF for genes of interest manually. Not entirely ideal but for a rough summary
# plot but this is going to be sufficient. Can revisit once we've unified AF information across callers:

# For AFs, use VarDict if found, otherwise Mutect2
# rule choose_af_vcf:
#     input:
#         vcfs = lambda wc: expand(join(run.date_dir, f'{batch_by_name[wc.batch].name}-{{caller}}-annotated.vcf.gz'), caller=list(batch_by_name.values())[0].tumor.variantcallers)
#     output:
#         '{batch}/somatic/af/variants.vcf.gz'
#     run:
#         af_vcf = next((vcf for vcf in input.vcfs if vcf.endswith('-ff-annotated.vcf.gz')), None)
#         if not af_vcf:
#             print('Could not find VarDict VCF, falling back to Mutect2')
#             af_vcf = next((vcf for vcf in input.vcfs if vcf.endswith('-mutect2-annotated.vcf.gz')), None)
#             assert af_vcf, 'Could not find Mutect2 VCF too'
#         shell('bcftools view -f.,PASS ' + af_vcf + ' -Oz -o {output}')

# Subset to GiaB confident intervals
rule subset_to_giab:
    input:
        vcf = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS.vcf.gz',
    params:
        regions = hpc.get_ref_file(run.genome_build, key=['hmf_giab_conf'])
    output:
        'work/{batch}/rmd/afs/' + run.somatic_caller + '-confident.vcf.gz'
    group: "rmd"
    shell:
        'bcftools view {input.vcf} -R {params.regions} -Oz -o {output}'

# Split multiallelics to avoid R parsing issues
rule split_multiallelic:
    input:
        vcf = 'work/{batch}/rmd/afs/' + run.somatic_caller + '-confident.vcf.gz',
        ref_fa = hpc.get_ref_file(run.genome_build, key='fa')
    output:
        'work/{batch}/rmd/afs/' + run.somatic_caller + '-confident-singleallelic.vcf.gz'
    group: "rmd"
    shell:
        'bcftools annotate -x ^INFO/TUMOR_AF {input.vcf} -Ob | '
        'bcftools norm -m \'-\' -Oz -f {input.ref_fa} -o {output} && tabix -p vcf {output}'

rule afs:
    input:
        'work/{batch}/rmd/afs/' + run.somatic_caller + '-confident-singleallelic.vcf.gz'
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name
    output:
        'work/{batch}/rmd/afs/af_tumor.txt'
    group: "rmd"
    shell:
        'bcftools view {input} -s {params.tumor_name} -Ou | '
        '(printf "af\n"; bcftools query -f "%INFO/TUMOR_AF\\n") > {output} && test -e {output}'

# Intersect with cancer key genes CDS for a table in Rmd
rule afs_keygenes:
    input:
        vcf = 'work/{batch}/rmd/afs/' + run.somatic_caller + '-confident-singleallelic.vcf.gz',
        bed = get_key_genes_bed(run.genome_build, coding_only=True),
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name
    output:
        'work/{batch}/rmd/afs/af_tumor_keygenes.txt'
    group: "rmd"
    shell:
        'bcftools view -f .,PASS {input.vcf} -s {params.tumor_name} -Ov'
        ' | bedtools intersect -a stdin -b {input.bed} -header'
        ' | (printf "chrom\\tpos\\tid\\tref\\talt\\taf\\n" ; bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/TUMOR_AF\\n")'
        ' > {output} && test -e {output}'

## Mutational signatures VCF
#
# Somatic calls to be submitted to CGI. Sharing functionality does not seem to work at this point.
# Somatic calls only include variants that `PASS` so no additional filtering required.
#
# Finally, for the local analysis with MutationalPatterns generate UCSC-versions (hg19) of the somatic calls:
rule somatic_to_hg19:
    input:
        vcf = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS.vcf.gz',
    output:
        'work/{batch}/rmd/' + run.somatic_caller + '-with_chr_prefix.vcf'
    group: "rmd"
    run:
        if run.genome_build == 'GRCh37':
            shell('gunzip -c {input}'
            " | py -x \"x.replace('##contig=<ID=', '##contig=<ID=chr') if x.startswith('#') else 'chr' + x\""
            " | py -x \"x.replace('chrMT', 'chrM')\""
            ' | grep -v chrG > {output}')
        else:
            shell('gunzip -c {input} > {output}')


# select chr, start, end, gene, min/max_cn, TranscriptID and ChromosomeBand
#rule rmd_purple_cnv:
#    input:
#        purple_cnv = rules.purple_run.output.gene_cnv,
#    output:
#        'work/{batch}/rmd/purple.tsv'
#    group: "rmd"
#    shell:
#        'cut -f1-6,11,13 {input} > {output}'

## Running Rmarkdown
rule cancer_report:
    input:
        rmd_files_dir        = join(package_path(), 'rmd_files'),
        key_genes            = get_key_genes(),
        af_global            = rules.afs.output[0],
        af_keygenes          = rules.afs_keygenes.output[0],
        somatic_snv          = rules.somatic_to_hg19.output[0],
        somatic_sv           = rules.prep_sv_tsv.output[0],
        purple_gene_cnv      = rules.purple_run.output.gene_cnv,
        purple_cnv           = rules.purple_run.output.cnv,
        purple_purity        = rules.purple_run.output.purity,
        purple_qc            = rules.purple_run.output.qc,

        purple_circos_png    = rules.purple_run.output.circos_png,
        purple_input_png     = rules.purple_run.output.input_png,
        purple_cn_png        = rules.purple_run.output.cn_png,
        purple_ma_png        = rules.purple_run.output.ma_png,
        purple_purity_png    = rules.purple_run.output.purity_png,
        purple_segment_png   = rules.purple_run.output.segment_png,
        purple_clonality_png = rules.purple_run.output.clonality_png,
        purple_ploidy_png    = rules.purple_run.output.ploidy_png,
        purple_rainfall_png  = rules.purple_run.output.rainfall_png,
        purple_baf_png       = rules.purple_circos_baf.output.png,

        purple_version       = rules.purple_run.output.version_build,
        amber_version        = rules.purple_amber.output.version_build,
        cobalt_version       = rules.purple_cobalt.output.version_build,

    params:
        report_rmd = 'cancer_report.Rmd',
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
        work_dir = os.getcwd(),
        output_file = lambda wc, output: join(os.getcwd(), output[0]),
        rmd_genome_build = 'hg19' if run.genome_build in ['GRCh37', 'hg19'] else run.genome_build,
        af_global           = lambda wc, input: abspath(input.af_global),
        af_keygenes         = lambda wc, input: abspath(input.af_keygenes),
        somatic_snv         = lambda wc, input: abspath(input.somatic_snv),
        somatic_sv          = lambda wc, input: abspath(input.somatic_sv),
        purple_gene_cnv     = lambda wc, input: abspath(input.purple_gene_cnv),
        purple_cnv          = lambda wc, input: abspath(input.purple_cnv),
        purple_purity       = lambda wc, input: abspath(input.purple_purity),
        purple_qc           = lambda wc, input: abspath(input.purple_qc),
        purple_version      = lambda wc, input: abspath(input.purple_version),
        amber_version       = lambda wc, input: abspath(input.amber_version),
        cobalt_version      = lambda wc, input: abspath(input.cobalt_version),
    output:
        report_html = '{batch}/{batch}_cancer_report.html',
        rmd_tmp_dir = directory('work/{batch}/rmd/rmd_files'),
    group: "rmd"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000
        # TODO: memory based on the mutation number. E.g. over 455k tumor mutations need over 5G
    run:
        shell('cp -r {input.rmd_files_dir} {output.rmd_tmp_dir}')
        shell('mkdir -p {output.rmd_tmp_dir}/img')
        for img_path in [
            input.purple_circos_png,
            input.purple_input_png,
            input.purple_cn_png,
            input.purple_ma_png,
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
purple_version='{params.purple_version}', \
amber_version='{params.amber_version}', \
cobalt_version='{params.cobalt_version}' \
))" ; \
cd {params.work_dir} ; \
""")


rule purple_bcbio_stats:
    input:
        rmd = verify_file(join(package_path(), 'rmd_files', 'purple_bcbio_umccrise.Rmd'), is_critical=True),
        purple_umccrise_files = expand('work/{batch}/purple/{batch}.purple.gene.cnv', batch=batch_by_name.keys()),
        bcbio_workdir = run.work_dir,
        key_genes = get_key_genes(),
    output:
        'purple_stats.html'
    params:
        rmd_tmp = 'work/purple/purple_bcbio_umccrise.Rmd',
        work_dir = abspath('work/purple'),
        output_file = lambda wc, output: abspath(output[0]),
    run:
        purple_bcbio_files = glob.glob(join(input.bcbio_workdir, 'structural/*/purple/purple/*.purple.gene.cnv'))
        assert(purple_bcbio_files)

        key_by_tumor_name = {
            b.tumor.name: bk for bk, b in batch_by_name.items()
        }
        safe_mkdir(params.workdir)
        for fn in purple_bcbio_files:
            print(fn)
            tumor_name = basename(fn).split('.purple.gene.cnv')[0]
            if tumor_name not in key_by_tumor_name:
                continue
            key = key_by_tumor_name[tumor_name]
            shell('cut -f1-5 ' + fn + ' > ' + join(params.workdir, f'bcbio_{key}.purple.gene.cnv'))
            print('Copying bcbio purple file to :', join(params.workdir, f'bcbio_{key}.purple.gene.cnv'))
        for fn in input.purple_umccrise_files:
            shell('cut -f1-5 ' + fn + ' > ' + join(params.workdir, f'umccrise_' + basename(fn)))
        shell(conda_cmd.format('cancer_report') + \
"""
cp {input.rmd} {params.rmd_tmp} && \
Rscript -e "rmarkdown::render('{params.rmd_tmp}',\
output_file='{params.output_file}', \
params=list( \
key_genes='{input.key_genes}', \
workdir='{params.work_dir}' \
))"
""")


#############

rule rmd:
    input:
        expand(rules.cancer_report.output[0], batch=batch_by_name.keys())
    output:
        temp(touch('log/rmd.done'))

