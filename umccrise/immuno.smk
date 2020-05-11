## Cancer gene coverage

from os.path import dirname

localrules: immuno


rule extract_hla_contigs:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
    output:
        'work/{batch}/immuno/{phenotype}_hla_contigs.txt'
    shell:
        "samtools view -H {input.bam} | "
        "grep 'SN:HLA' | cut -f2 | sed 's/SN://' > {output}"


rule extact_hla_bam:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        hla_contigs = rules.extract_hla_contigs.output[0],
    output:
        bam = 'work/{batch}/immuno/{phenotype}_hla.bam'
    params:
        hla_region = {'hg38': 'chr6:29942470-31357179',
                      'GRCh37': '6:29910247-31324989'}[run.genome_build]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    threads:
        threads_per_sample
    group: 'hla'
    run:
        with open(input.hla_contigs) as f:
            contigs = f.readlines()
        if contigs:
            # The reference with HLA as separate contigs was used. Subsetting to HLA contigs:
            shell('samtools view {input.bam} $(cat {input.hla_contigs}) -Obam > {output.bam}')
        else:
            # The refernece without HLA contigs was used. Slicing chromosome 6:
            shell('sambamba slice {input.bam} {params.hla_region} > {output.bam}')


rule hla_fastq:
    input:
        bam = rules.extact_hla_bam.output.bam
    output:
        fastq = 'work/{batch}/immuno/{phenotype}_hla.fastq'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    threads:
        threads_per_sample
    group: 'hla'
    shell:
        'samtools fastq -@{threads} {input.bam} -o {output.fastq}'


rule run_optitype:
    input:
        fastq = rules.hla_fastq.output.fastq
    output:
        tsv = '{batch}/hla/{batch}_{phenotype}_result.tsv',
        pdf = '{batch}/hla/{batch}_{phenotype}_coverage_plot.pdf',
    params:
        out_dir = 'work/{batch}/hla/',
        prefix = '{batch}_{phenotype}',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000
    threads:
        threads_per_sample
    group: 'hla'
    benchmark:
        'benchmarks/{batch}/immuno/{batch}-{phenotype}-optitype.tsv'
    shell:
        'OptiTypePipeline.py -i {input.fastq} --dna '
        '--outdir {params.out_dir} --prefix {params.prefix}'


rule immuno:
    input:
        expand(rules.run_optitype.output, batch=batch_by_name.keys(), phenotype=['tumor', 'normal'])
    output:
        temp(touch('log/immuno.done'))














