from umccrise.patient_analysis import get_sig_rmd_file


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
# All done; file needs to be submitted to CGI and PCGR manually for now. 


rule restrict_to_confident:
    input:
        rules.cgi.output[0]
    params:
        regions = join(loc.hsapiens, loc.truth_sets['giab'][run.genome_build]['bed'])
    output:
        'umccrised/{batch}/somatic/rstudio/ensemble-ucsf-confident.vcf'
    shell:
        'bcftools view {input} -T {params.regions} -Oz -o {output}'


rule sig_rmd:
    input:
        af_freqs = rules.af_freqs.output[0],
        af_freqs_az300 = rules.af_freqs_az300.output[0],
        vcf = rules.restrict_to_confident.output[0]
    params:
        sig_rmd = get_sig_rmd_file(),
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
        workdir = os.getcwd(),
        output_file = join(os.getcwd(), 'umccrised/{batch}/somatic/rstudio/sig.html')
    output:
        'umccrised/{batch}/somatic/rstudio/sig.html'
    run:
        os.system("which Rscript")
        shell('Rscript -e "rmarkdown::render(\'{params.sig_rmd}\', '
        'output_file=\'{params.output_file}\', '
        'params=list('
        'af_freqs=\'{input.af_freqs}\', '
        'af_freqs_az300=\'{input.af_freqs_az300}\', '
        'vcf=\'{input.vcf}\', '
        'tumor_name=\'{params.tumor_name}\', '
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
    
