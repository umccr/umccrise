## Cancer gene coverage
from os.path import join, dirname

from ngs_utils.file_utils import which

localrules: immuno


"""
Using OptiType to predict HLA types. As input, extracting 3 categories of reads:
- Reads aligned to HLA hg38 contigs (however they barely contribute anything)
- Reads mapped to the chr6:29942470-31357179 region (significantly improves sensitivity)
- Unmapped reads (slightly improves sensitivity)

Also, adding razers2 step to filter reads to those taht align the HLA reference provided by OptiType.
This step is actually done by OptiType itself, however doing it separately slightly improves
sensitivity and also gives us more control over the pipeline flow and resources. 

Comparison of taking different sources of input for OptiType:

Dragen+umccrise, using only reads mapping to hg38 HLA choromosomes:
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*01:01 A*02:01 B*08:01 B*08:01                 22.0    21.802
N: 0       A*01:01 A*01:01 B*08:01 B*35:01 C*07:01 C*07:01 10.0    9.91

Dragen+umccrise, plus the reads mapped to the chr6 HLA region:
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 1435.0  1383.3400000000008
N: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 765.0   737.4600000000011

Dragen+umccrise, also with unmapped reads:
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 1681.0  1620.484000000002
N: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 929.0   895.5560000000014

Dragen+umccrise, with an extra intermediate razers2 step:
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 2384.0  2298.176
N: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:17 1302.0  1255.1180000000015

Dragen+umccrise, Running razers2 on read1 and read2 separately, separately fed to OptiType with `-i r1.fq r2.fq`:
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 703.0   677.6920000000002
N: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 373.0   359.57199999999983

bcbio standalone (internally it merges read1 and read1 together in one file before passing to OptiType):
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*02:01 A*01:01 B*07:02 B*08:01 C*07:02 C*07:01 1758.0  1710.5339999999999        
N: 0       A*02:01 A*01:01 B*07:02 B*08:01 C*07:01 C*07:17 1008.0  980.7739999999995

bcbio+umccrise, read1 and read2 razers'ed separately, separately fed to OptiType with `-i r1.fq r2.fq`:
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 685.0   660.3400000000001
N: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 359.0   346.0759999999999

bcbio+umccrise, umccrise for read extraction, but bcbio OptiType command. Same result as umccrise on bcbio with `-i r1.fq r2.fq`:
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*02:01 A*01:01 B*07:02 B*08:01 C*07:02 C*07:01 703.0   677.692

bcbio+umccrise, read1 and read2 merged and fed into optitype together:
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:02 2365.0  2279.8599999999997
N: 0       A*01:01 A*02:01 B*07:02 B*08:01 C*07:01 C*07:17 1263.0  1217.522000000001

Feeding read1 and read2 together or separately to OptiType makes difference. Just the matter of
OptiTypePipeline.py -v --dna -o unmerged -i tumor.1.fastq tumor.2.fastq -c config.ini
or
cat tumor.1.fastq tumor.2.fastq > tumor.fastq
OptiTypePipeline.py -v --dna -o unmerged -i tumor.fastq -c config.ini

What if we try -i -i ?
OptiTypePipeline.py -v --dna -o unmerged_ii -i tumor.1.fastq -i tumor.2.fastq -c config.ini
           A1      A2      B1      B2      C1      C2      Reads   Objective
T: 0       A*02:01 A*01:01 B*07:02 B*08:01 C*07:02 C*07:01 1105.0  1065.2200000000014

Different result again. Sticking to pre-merging fastq.
"""

HLA_READ_SOURCES = [
    'from_contigs',
    'from_chr6',
    # 'from_unmapped',
]


rule extract_hla_from_unmapped:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
    output:
        bam = 'work/{batch}/hla/from_unmapped/{phenotype}.bam'
    threads:
        threads_per_sample,
    resources:
        mem_mb=10000
    shell:
         "sambamba view -t{threads} -fbam -F 'unmapped' {input.bam} -o {output.bam}"


rule extract_hla_from_chr6:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
    output:
        bam = 'work/{batch}/hla/from_chr6/{phenotype}.bam'
    params:
        hla_region = {'hg38': 'chr6:29942470-31357179',
                      'GRCh37': '6:29910247-31324989'}[run.genome_build]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    group: 'hla'
    shell:
        'sambamba slice {input.bam} {params.hla_region} > {output.bam}'


rule extract_hla_contigs:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
    output:
        'work/{batch}/hla/from_contigs/{phenotype}_hla_contigs.txt'
    shell:
        "samtools view -H {input.bam} | "
        "grep 'SN:HLA' | cut -f2 | sed 's/SN://' > {output}"

rule extract_hla_from_hla_contigs:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        hla_contigs = rules.extract_hla_contigs.output[0],
    output:
        bam = 'work/{batch}/hla/from_contigs/{phenotype}.bam'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    group: 'hla'
    shell:
        'samtools view {input.bam} $(cat {input.hla_contigs}) -Obam > {output.bam}'

rule hla_fastq:
    input:
        bam = 'work/{batch}/hla/{read_source}/{phenotype}.bam'
    output:
        fastq1 = 'work/{batch}/hla/{read_source}/{phenotype}.1.fastq',
        fastq2 = 'work/{batch}/hla/{read_source}/{phenotype}.2.fastq',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    threads:
        threads_per_sample
    group: 'hla'
    shell:
        'samtools fastq -@{threads} {input.bam} -1 {output.fastq1} -2 {output.fastq2} -N'

rule hla_razers:
    input:
        fastq = 'work/{batch}/hla/{read_source}/{phenotype}.{num}.fastq',
        hla_ref = join(dirname(which('OptiTypePipeline.py')), 'data', 'hla_reference_dna.fasta')
    output:
        bam = 'work/{batch}/hla/razers_{read_source}/{phenotype}.{num}.bam'
    group: 'hla'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000
    threads:
        threads_per_sample
    shell:
        'razers3 {input.hla_ref} {input.fastq}'
        ' --percent-identity 95 --max-hits 1 --distance-range 0'
        ' --thread-count {threads} --output {output.bam}'

rule hla_razers_to_fastq:
    input:
        bam = 'work/{batch}/hla/razers_{read_source}/{phenotype}.{num}.bam',
    output:
        fastq = 'work/{batch}/hla/razers_{read_source}/{phenotype}.{num}.fastq'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    threads:
        threads_per_sample
    group: 'hla'
    shell:
        'samtools fastq -@{threads} {input.bam} > {output.fastq}'

rule merge_razers_fastqs:
    input:
        fastqs = expand('work/{{batch}}/hla/razers_{read_source}/{{phenotype}}.{num}.fastq',
                        read_source=HLA_READ_SOURCES, num=[1, 2]),
    output:
        fastq = 'work/{batch}/hla/optitype_input/{phenotype}.fastq'
    group: 'hla'
    shell:
        'cat {input} > {output.fastq}'

rule run_optitype:
    input:
        fastq1 = 'work/{batch}/hla/optitype_input/{phenotype}.fastq',
        # fastq2 = 'work/{batch}/hla/optitype_input/{phenotype}.2.fastq',
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
        'OptiTypePipeline.py -i {input.fastq1} --dna --verbose '
        '--outdir {params.out_dir} --prefix {params.prefix}'


rule immuno:
    input:
        expand(rules.run_optitype.output, batch=batch_by_name.keys(), phenotype=['tumor', 'normal'])
    output:
        temp(touch('log/immuno.done'))








