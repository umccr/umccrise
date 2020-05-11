## Cancer gene coverage

from os.path import dirname

localrules: immuno


rule extract_hla_contigs:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
    output:
        'work/{batch}/hla/{phenotype}_hla_contigs.txt'
    shell:
        "samtools view -H {input.bam} | "
        "grep 'SN:HLA' | cut -f2 | sed 's/SN://' > {output}"


rule extract_hla_from_chr6:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
    output:
        bam = 'work/{batch}/hla/from_contigs/{phenotype}_hla.bam'
    params:
        hla_region = {'hg38': 'chr6:29942470-31357179',
                      'GRCh37': '6:29910247-31324989'}[run.genome_build]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    threads:
        threads_per_sample
    group: 'hla'
    shell:
        'sambamba slice {input.bam} {params.hla_region} > {output.bam}'


rule extract_hla_from_hla_contigs:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        hla_contigs = rules.extract_hla_contigs.output[0],
    output:
        bam = 'work/{batch}/hla/from_chr6/{phenotype}_hla.bam'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    threads:
        threads_per_sample
    group: 'hla'
    shell:
        'samtools view {input.bam} $(cat {input.hla_contigs}) -Obam > {output.bam}'


rule hla_fastq_from_hla_contigs:
    input:
        bam = rules.extract_hla_from_hla_contigs.output.bam,
    output:
        fastq = 'work/{batch}/hla/from_contigs/{phenotype}_hla.fastq'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    threads:
        threads_per_sample
    group: 'hla'
    shell:
        'samtools fastq -@{threads} {input.bam} -o {output.fastq}'

rule hla_fastq_from_chr6:
    input:
        bam = rules.extract_hla_from_chr6.output.bam,
    output:
        fastq = 'work/{batch}/hla/from_chr6/{phenotype}_hla.fastq'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    threads:
        threads_per_sample
    group: 'hla'
    shell:
        'samtools fastq -@{threads} {input.bam} -o {output.fastq}'


rule merge_hla_fastq:
    input:
        rules.hla_fastq_from_hla_contigs.output.fastq,
        rules.hla_fastq_from_chr6.output.fastq,
    output:
        fastq = 'work/{batch}/hla/{phenotype}_hla.fastq'
    group: 'hla'
    shell:
        'cat {input} > {output.fastq}'


rule run_optitype:
    input:
        fastq = rules.merge_hla_fastq.output.fastq
    output:
        tsv = '{batch}/hla/{batch}_{phenotype}_result.tsv',
        pdf = '{batch}/hla/{batch}_{phenotype}_coverage_plot.pdf',
    params:
        out_dir = '{batch}/hla',
        prefix = '{batch}_{phenotype}',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000
    threads:
        threads_per_sample
    group: 'hla'
    benchmark:
        'benchmarks/{batch}/hla/{batch}-{phenotype}-optitype.tsv'
    shell:
        'OptiTypePipeline.py -i {input.fastq} --dna '
        '--outdir {params.out_dir} --prefix {params.prefix}'


rule immuno:
    input:
        expand(rules.run_optitype.output, batch=batch_by_name.keys(), phenotype=['tumor', 'normal'])
    output:
        temp(touch('log/immuno.done'))














