from umccrise.patient_analysis import get_sig_rmd_file, get_signatures_probabilities


## Allelic frequencies
# AF is not yet integrated into PCGR or CGI. We can extract those from VarDict or Mutect2 for plotting purposes,
# but still need to look up AF for genes of interest manually. Not entirely ideal but for a rough summary
# plot but this is going to be sufficient. Can revisit once we've unified AF information across callers:

# For AFs, use VarDict if found, otherwise Mutect2
rule choose_vcf:
    input:
        vcfs = lambda wc: expand(join(run.date_dir, f'{batch_by_name[wc.batch].name}-{{caller}}-annotated.vcf.gz'), caller=list(batch_by_name.values())[0].tumor.variantcallers)
    output:
        'umccrised/{batch}/somatic/af/variants.vcf.gz'
    run:
        af_vcf = next((vcf for vcf in input.vcfs if vcf.endswith('-ff-annotated.vcf.gz')), None)
        if not af_vcf:
            print('Could not find VarDict VCF, falling back to Mutect2')
            af_vcf = next((vcf for vcf in input.vcfs if vcf.endswith('-mutect2-annotated.vcf.gz')), None)
            assert af_vcf, 'Could not find Mutect2 VCF too'
        shell('cp ' + af_vcf + ' {output}')

# Subset to GiaB confident intervals
rule restrict_to_confident:
    input:
        vcf = rules.choose_vcf.output[0]
    params:
        regions = join(loc.hsapiens, loc.truth_sets['giab'][run.genome_build]['bed'])
    output:
        'umccrised/{batch}/somatic/af/variants.confident.vcf.gz'
    shell:
        'bcftools view {input.vcf} -T {params.regions} -Oz -o {output}'

# Split multiallelics to avoid R parsing issues
rule normalise_for_afs:
    input:
        vcf = rules.restrict_to_confident.output[0],
        ref_fa = ref_fa
    output:
        'umccrised/{batch}/somatic/af/variants.confident.normalised.vcf.gz'
    shell:
        "bcftools annotate -x ^INFO/ANN -Ob {input.vcf} | bcftools norm -m '-' -Oz -f {input.ref_fa} -o {output} && tabix -p vcf {output}"

rule afs:
    input:
        rules.normalise_for_afs.output[0]
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name
    output:
        'umccrised/{batch}/somatic/af/af_tumor.txt'
    shell:
        'bcftools view -f .,PASS {input} -s {params.tumor_name} -Ou | bcftools query -f "[%AF]\\n" > {output}'

rule afs_az300:
    input:
        vcf = rules.normalise_for_afs.output[0],
        az300 = az300
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name
    output:
        'umccrised/{batch}/somatic/af/af_tumor_az300.txt'
    shell:
        'bcftools view -f .,PASS {input.vcf} -s {params.tumor_name} -Ov'
        ' | bedtools intersect -a stdin -b {input.az300} -header'
        ' | bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%AF]\\t%INFO/ANN\n" > {output}'

## Mutational signatures VCF
# Ensemble calls to be submitted to CGI. Sharing functionality does not seem to work at this point. Ensemble calls only include variants that `PASS` so no additional filtering required.
# Finally, for the local analysis with MutationalPatterns generate UCSC-versions of the ensemble calls:
rule cgi:
    input:
        'umccrised/{batch}/somatic/pcgr/ensemble.vcf.gz'
    output:
        'umccrised/{batch}/somatic/rstudio/ensemble-ucsf.vcf'
    shell:
        'gunzip -c {input}'
        " | py -x \"x.replace('##contig=<ID=', '##contig=<ID=chr') if x.startswith('#') else 'chr' + x\""
        " | py -x \"x.replace('chrMT', 'chrM')\""
        ' | grep -v chrG > {output}'


## Running Rmarkdown
rule sig_rmd:
    input:
        afs = rules.afs.output[0],
        afs_az300 = rules.afs_az300.output[0],
        vcf = rules.cgi.output[0],
        sig_rmd = get_sig_rmd_file(),
        sig_probs = get_signatures_probabilities()
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
        workdir = os.getcwd(),
        output_file = join(os.getcwd(), 'umccrised/{batch}/somatic/rstudio/sig.html')
    output:
        'umccrised/{batch}/somatic/rstudio/sig.html'
    run:
        os.system("which Rscript")
        shell('Rscript -e "rmarkdown::render(\'{input.sig_rmd}\', '
        'output_file=\'{params.output_file}\', '
        'params=list('
        'af_freqs=\'{input.afs}\', '
        'af_freqs_az300=\'{input.afs_az300}\', '
        'vcf=\'{input.vcf}\', '
        'tumor_name=\'{params.tumor_name}\', '
        'sig_probs=\'{input.sig_probs}\', '
        'workdir=\'{params.workdir}\''
        '))"'
        ' && if [ -e x.txt ] ; then exit 1; fi')


rule sig:
    input:
        expand(rules.sig_rmd.output[0], batch=batch_by_name.keys())
    output:
        'umccrised/.snakemake/sig.done'
    shell:
        'touch {output}'
    
