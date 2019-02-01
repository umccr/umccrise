from os.path import join

from ngs_utils.logger import warn
from ngs_utils.reference_data import get_key_genes, get_key_genes_bed
from ngs_utils.file_utils import safe_mkdir
from umccrise import get_sig_rmd_file, get_signatures_probabilities
import glob


# localrules: purple_bcbio_stats, rmd
localrules: split_multiallelic, subset_to_giab, afs, afs_keygenes, somatic_to_hg19, rmd_purple, purple_bcbio_stats, rmd


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
        vcf = rules.somatic_vcf_filter_pass.output.vcf
    params:
        regions = truth_regions
    output:
        'work/{batch}/rmd/afs/' + run.somatic_caller + '-confident.vcf.gz'
    # group: "rmd"
    shell:
        'bcftools view {input.vcf} -T {params.regions} -Oz -o {output}'

# Split multiallelics to avoid R parsing issues
rule split_multiallelic:
    input:
        vcf = 'work/{batch}/rmd/afs/' + run.somatic_caller + '-confident.vcf.gz',
        ref_fa = ref_fa
    output:
        'work/{batch}/rmd/afs/' + run.somatic_caller + '-confident-singleallelic.vcf.gz'
    # group: "rmd"
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
    # group: "rmd"
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
    # group: "rmd"
    shell:
        'bcftools view -f .,PASS {input.vcf} -s {params.tumor_name} -Ov'
        ' | bedtools intersect -a stdin -b {input.bed} -header'
        ' | (printf "chrom\tpos\tid\tref\talt\taf\n" ; bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/TUMOR_AF\\n")'
        ' > {output} && test -e {output}'

## Mutational signatures VCF
#
# Somatic calls to be submitted to CGI. Sharing functionality does not seem to work at this point.
# Somatic calls only include variants that `PASS` so no additional filtering required.
#
# Finally, for the local analysis with MutationalPatterns generate UCSC-versions (hg19) of the somatic calls:
rule somatic_to_hg19:
    input:
        rules.somatic_vcf_filter_pass.output.vcf
    output:
        'work/{batch}/rmd/' + run.somatic_caller + '-with_chr_prefix.vcf'
    # group: "rmd"
    run:
        if run.genome_build == 'GRCh37':
            shell('gunzip -c {input}'
            " | py -x \"x.replace('##contig=<ID=', '##contig=<ID=chr') if x.startswith('#') else 'chr' + x\""
            " | py -x \"x.replace('chrMT', 'chrM')\""
            ' | grep -v chrG > {output}')
        else:
            shell('gunzip -c {input} > {output}')


rule rmd_purple_cnv:
    input:
        'work/{batch}/purple/{batch}.purple.gene.cnv',
    output:
        'work/{batch}/rmd/purple.tsv'
    # group: "rmd"
    shell:
        'cut -f1-6 {input} > {output}'


# ## Running Rmarkdown
# rule sig_rmd:
#     input:
#         afs = rules.afs.output[0],
#         afs_keygenes = rules.afs_keygenes.output[0],
#         vcf = rules.somatic_to_hg19.output[0],
#         sv = rules.prep_sv_tsv.output[0],
#         sig_rmd = get_sig_rmd_file(),
#         sig_probs = get_signatures_probabilities(),
#         key_genes = get_key_genes(),
#         manta_vcf = rules.filter_sv_vcf.output[0],
#         purple = rules.rmd_purple.output[0],
#     params:
#         rmd_tmp = 'work/{batch}/rmd/sig.Rmd',
#         tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
#         workdir = os.getcwd(),
#         output_file = lambda wc, output: join(os.getcwd(), output[0]),
#         rmd_genome_build = 'hg19' if run.genome_build in ['GRCh37', 'hg19'] else run.genome_build
#     output:
#         '{batch}/{batch}-rmd_report.html'
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 10000
#         # TODO: memory based on the mutation number. E.g. over 455k tumor mutations need over 5G
#     shell: """
# cp {input.sig_rmd} {params.rmd_tmp} &&
# Rscript -e "rmarkdown::render('{params.rmd_tmp}', \
# output_file='{params.output_file}', \
# params=list( \
# af_freqs='{input.afs}', \
# af_freqs_keygenes='{input.afs_keygenes}', \
# vcf_fname='{input.vcf}', \
# sv_fname='{input.sv}', \
# manta_vcf='{input.manta_vcf}', \
# tumor_name='{params.tumor_name}', \
# sig_probs='{input.sig_probs}', \
# key_genes='{input.key_genes}', \
# purple='{input.purple}', \
# workdir='{params.workdir}', \
# genome_build='{params.rmd_genome_build}' \
# ))"
# """

## Running Rmarkdown
rule bookdown_report:
    input:
        rmd_files_dir = join(package_path(), 'rmd_files'),
        sig_probs = join(package_path(), 'rmd_files', 'signatures_probabilities.txt'),
        key_genes = get_key_genes(),
        afs                 = rules.afs.output[0],
        afs_keygenes        = rules.afs_keygenes.output[0],
        vcf                 = rules.somatic_to_hg19.output[0],
        sv                  = rules.prep_sv_tsv.output[0],
        manta_vcf           = rules.filter_sv_vcf.output[0],
        purple_gene_cnv     = rules.rmd_purple_cnv.output[0],
        purple_cnv          = rules.purple_run.output.cnv,
        purple_germline_cnv = rules.purple_run.output.germline_cnv,
        purple_purity       = rules.purple_run.output.purity,
        purple_circos_png   = rules.purple_run.output.circos_png,
        purple_input_png    = rules.purple_run.output.input_png,
        purple_cn_png       = rules.purple_run.output.cn_png,
        purple_ma_png       = rules.purple_run.output.ma_png,
        purple_variant_png  = rules.purple_run.output.variant_png,
        purple_baf_png      = rules.purple_circos_baf.output.png,
    params:
        index_rmd = 'index.Rmd',
        # bookdown_yml = 'work/{batch}/rmd/rmd_files/_bookdown.yml',
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
        work_dir = os.getcwd(),
        output_file = lambda wc, output: join(os.getcwd(), output[0]),
        rmd_genome_build = 'hg19' if run.genome_build in ['GRCh37', 'hg19'] else run.genome_build,
        afs                 = lambda wc, input: abspath(input.afs),
        afs_keygenes        = lambda wc, input: abspath(input.afs_keygenes),
        vcf                 = lambda wc, input: abspath(input.vcf),
        sv                  = lambda wc, input: abspath(input.sv),
        manta_vcf           = lambda wc, input: abspath(input.manta_vcf),
        purple_gene_cnv     = lambda wc, input: abspath(input.purple_gene_cnv),
        purple_cnv          = lambda wc, input: abspath(input.purple_cnv),
        purple_germline_cnv = lambda wc, input: abspath(input.purple_germline_cnv),
        purple_purity       = lambda wc, input: abspath(input.purple_purity),
    output:
        report_html = '{batch}/{batch}_book.html',
        rmd_tmp_dir = directory('work/{batch}/rmd/rmd_files'),
    # group: "rmd"
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
            input.purple_variant_png,
            input.purple_baf_png,
        ]:
            shell('cp ' + img_path + ' {output.rmd_tmp_dir}/img/')
        shell("""
cd {output.rmd_tmp_dir} && \
Rscript -e "library(bookdown); bookdown::render_book('{params.index_rmd}', \
params=list( \
tumor_name='{params.tumor_name}', \
batch_name='{wildcards.batch}', \
genome_build='{params.rmd_genome_build}', \
sig_probs='{input.sig_probs}', \
key_genes='{input.key_genes}', \
af_freqs='{params.afs}', \
af_freqs_keygenes='{params.afs_keygenes}', \
vcf_fname='{params.vcf}', \
sv_fname='{params.sv}', \
manta_vcf='{params.manta_vcf}', \
purple_gene_cnv='{params.purple_gene_cnv}', \
purple_cnv='{params.purple_cnv}', \
purple_germline_cnv='{params.purple_germline_cnv}', \
purple_purity='{params.purple_purity}' \
))" ; \
cd {params.work_dir} ; \
""")
        shell('mv {output.rmd_tmp_dir}/_book_umccrised/_main.html {output.report_html}')

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
        shell("""
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
        expand(rules.bookdown_report.output[0], batch=batch_by_name.keys())
    output:
        temp(touch('log/rmd.done'))

