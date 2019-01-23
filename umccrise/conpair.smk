## Cancer gene coverage

localrules: conpair


rule run_conpair:
    input:
        tumor_bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
    output:
        directory('{batch}/conpair/concordance'),
        directory('{batch}/conpair/contamination'),
    threads: 2
    params:
        genome = run.genome_build,
        out_dir = '{batch}/conpair',
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
        normal_name = lambda wc: batch_by_name[wc.batch].normal.name,
    shell:
        'conpair -T {input.tumor_bam} -N {input.normal_bam} -g {params.genome} -j {threads} '
        '-o {params.out_dir} -tn {params.tumor_name} -nn {params.normal_name}'


rule conpair:
    input:
        expand(rules.run_conpair.output, batch=batch_by_name.keys())
    output:
        temp(touch('log/conpair.done'))
