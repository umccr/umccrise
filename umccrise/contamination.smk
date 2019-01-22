## Cancer gene coverage

localrules: contamination


rule conpair:
    input:
        tumor_bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
    output:
        '{batch}/conpair/tumor_normal_concordance.txt',
        '{batch}/conpair/tumor_normal_contamination.txt',
    threads: 2
    params:
        genome = run.genome_build,
        out_dir = '{batch}/conpair'
    shell:
        'conpair -T {input.tumor_bam} -N {input.normal_bam} -g {params.genome} -j {threads} '
        '-o {params.out_dir} -tn tumor -nn normal'


rule contamination:
    input:
        expand(rules.conpair.output[0], batch=batch_by_name.keys())
    output:
        temp(touch('log/contamination.done'))
