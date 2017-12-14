
rule sig_rmd:
    input:
        af_freqs = rules.af_freqs.output[0],
        af_freqs_az300 = rules.af_freqs_az300.output[0],
        vcf = rules.cgi.output[0]
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
        expand(rules.sig_rmd.output[0], batch=batch_by_name.values())
    output:
        'umccrised/.snakemake/sig.done'
    shell:
        'touch {output}'
    
