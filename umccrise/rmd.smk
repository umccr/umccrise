from os.path import join
from ngs_utils.reference_data import get_suppressors
from umccrise import get_sig_rmd_file, get_signatures_probabilities


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
        vcf = rules.pcgr_somatic_vcf.output.vcf
    params:
        regions = truth_regions
    output:
        'work/{batch}/rmd/afs/ensemble-confident.vcf.gz'
    group: "subset_for_af"
    shell:
        'bcftools view {input.vcf} -T {params.regions} -Oz -o {output}'

# Split multiallelics to avoid R parsing issues
rule split_multiallelic:
    input:
        vcf = 'work/{batch}/rmd/afs/ensemble-confident.vcf.gz',
        ref_fa = ref_fa
    output:
        'work/{batch}/rmd/afs/ensemble-confident-singleallelic.vcf.gz'
    group: "subset_for_af"
    shell:
        "bcftools annotate -x ^INFO/ANN,^INFO/TUMOR_AF -Ob {input.vcf} | "
        "bcftools norm -m '-' -Oz -f {input.ref_fa} -o {output} && tabix -p vcf {output}"

rule afs:
    input:
        'work/{batch}/rmd/afs/ensemble-confident-singleallelic.vcf.gz'
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name
    output:
        'work/{batch}/rmd/afs/af_tumor.txt'
    group: "subset_for_af"
    shell:
        'bcftools view {input} -s {params.tumor_name} -Ou | '
        'bcftools query -f "%INFO/TUMOR_AF\\n" > {output} && test -e {output}'

rule afs_az300:
    input:
        vcf = 'work/{batch}/rmd/afs/ensemble-confident-singleallelic.vcf.gz',
        az300 = key_genes_bed
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name
    output:
        'work/{batch}/rmd/afs/af_tumor_az300.txt'
    group: "subset_for_af"
    shell:
        'bcftools view -f .,PASS {input.vcf} -s {params.tumor_name} -Ov'
        ' | bedtools intersect -a stdin -b {input.az300} -header'
        ' | bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/TUMOR_AF\\t%INFO/ANN\\n"'
        ' > {output} && test -e {output}'

## Mutational signatures VCF
#
# Ensemble calls to be submitted to CGI. Sharing functionality does not seem to work at this point.
# Ensemble calls only include variants that `PASS` so no additional filtering required.
#
# Finally, for the local analysis with MutationalPatterns generate UCSC-versions (hg19) of the ensemble calls:
rule somatic_to_hg19:
    input:
        rules.pcgr_somatic_vcf.output.vcf
    output:
        'work/{batch}/rmd/ensemble-with_chr_prefix.vcf'
    group: "subset_for_af"
    run:
        if run.genome_build == 'GRCh37':
            shell('gunzip -c {input}'
            " | py -x \"x.replace('##contig=<ID=', '##contig=<ID=chr') if x.startswith('#') else 'chr' + x\""
            " | py -x \"x.replace('chrMT', 'chrM')\""
            ' | grep -v chrG > {output}')
        else:
            shell('gunzip -c {input} > {output}')


## Running Rmarkdown
rule sig_rmd:
    input:
        afs = rules.afs.output[0],
        afs_az300 = rules.afs_az300.output[0],
        vcf = rules.somatic_to_hg19.output[0],
        sv = rules.prep_sv_tsv.output[0],
        sig_rmd = get_sig_rmd_file(),
        sig_probs = get_signatures_probabilities(),
        suppressors = get_suppressors(),
        cnvkit_calls = lambda wc: join(batch_by_name[wc.batch].tumor.dirpath, f'{batch_by_name[wc.batch].name}-cnvkit-call.cns'),
        manta_vcf = rules.filter_sv_vcf.output[0]
    params:
        rmd_tmp = 'work/{batch}/rmd/sig.Rmd',
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
        workdir = os.getcwd(),
        output_file = lambda wc, output: join(os.getcwd(), output[0]),
        rmd_genome_build = 'hg19' if run.genome_build in ['GRCh37', 'hg19'] else run.genome_build
    output:
        '{batch}/{batch}-rmd_report.html'
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    shell: """cp {input.sig_rmd} {params.rmd_tmp} && 
Rscript -e "rmarkdown::render('{params.rmd_tmp}', \
output_file='{params.output_file}', \
params=list( \
af_freqs='{input.afs}', \
af_freqs_az300='{input.afs_az300}', \
vcf_fname='{input.vcf}', \
sv_fname='{input.sv}', \
cnvkit_calls='{input.cnvkit_calls}', \
manta_vcf='{input.manta_vcf}', \
tumor_name='{params.tumor_name}', \
sig_probs='{input.sig_probs}', \
suppressors='{input.suppressors}', \
workdir='{params.workdir}', \
genome_build='{params.rmd_genome_build}' \
))"
"""


#############

rule rmd:
    input:
        expand(rules.sig_rmd.output[0], batch=batch_by_name.keys())
    output:
        temp(touch('log/rmd.done'))