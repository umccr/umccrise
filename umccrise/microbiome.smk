from collections import defaultdict
from os.path import join, basename
import yaml
from cyvcf2 import VCF
import csv
from ngs_utils.file_utils import verify_file
from ngs_utils.utils import update_dict
from reference_data import api as refdata


localrules: microbiome


"""
Working with unmapped reads

Tools:
- kraken
- blast
- mash scren
- rkmh https://github.com/edawson/rkmh rkmh performs identification of individual reads, identity-based
  read filtering, and alignment-free variant calling using MinHash (as implemented in Mash). It is compatible 
  with Mash and sourmash via JSON exchange.

"""

if 'microbiome' in stages:
    rule extract_unmapped:
        input:
            host_bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
        output:
            host_bam_namesorted = 'work/{batch}/microbiome/step1_unmapped.namesorted.bam',
        threads:
            threads_per_batch,
        resources:
            mem_mb=10000
        shell:
             "sambamba view -t{threads} -fbam -F 'unmapped' {input.host_bam}"
             " | samtools sort -n -@{threads} -Obam -o {output.host_bam_namesorted}"


    rule unmapped_to_fastq:
        input:
            host_bam_namesorted = rules.extract_unmapped.output.host_bam_namesorted,
        output:
            fq1 = 'work/{batch}/microbiome/step2_host_unmapped.R1.fq',
            fq2 = 'work/{batch}/microbiome/step2_host_unmapped.R2.fq',
            single = 'work/{batch}/microbiome/step2_host_unmapped.single.fq',
        threads:
            threads_per_batch,
        resources:
            mem_mb=4000
        shell:
            "samtools fastq -@{threads} {input.host_bam_namesorted} "
            "-1 {output.fq1} -2 {output.fq2} -s {output.single}"


    rule merge_reads:
        input:
            'work/{batch}/microbiome/step2_host_unmapped.R1.fq',
            'work/{batch}/microbiome/step2_host_unmapped.R2.fq',
            'work/{batch}/microbiome/step2_host_unmapped.single.fq',
        output:
            fq = 'work/{batch}/microbiome/step3_unmapped_merged.fq',
        shell:
            'cat {input} > {output.fq}'


    rule run_mash_screen:
        input:
            fq = rules.merge_reads.output.fq,
            refseq_msh = refdata.get_ref_file(key = 'refseq_msh'),
        output:
            screen_tab = 'work/{batch}/microbiome/step3_screen.tab',
        threads:
            threads_per_batch,
        resources:
            mem_mb=4000
        shell:
            conda_cmd.format('microbiome') +
            'mash screen {input.refseq_msh} {input.fq} | sort -gr > {output.screen_tab}'


    rule align_to_4bacteria:
        input:
            seq = '/g/data/gx8/projects/Saveliev_APGI/LTS/krakenuniq/db_4bacterias/library/combined.fasta',
            r1 = 'work/{batch}/microbiome/step2_host_unmapped.R1.fq',
            r2 = 'work/{batch}/microbiome/step2_host_unmapped.R2.fq',
        output:
            'work/{batch}/microbiome/step4_bwa_to_4bacteria.bam'
        shell:
            'bwa mem {input.seq} {input.r1} {input.r2} | samtools sort -Obam -o {output} '
            '&& samtools index {output}'


    rule microbiome:
        input:
            expand(rules.run_mash_screen.output, batch=batch_by_name.keys()),
            # expand(rules.align_to_4bacteria.output, batch=batch_by_name.keys()),
        output:
            temp(touch('log/microbiome.done'))













