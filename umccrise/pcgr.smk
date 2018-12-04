"""
PCGR
-------------
Prepare somatic, germline variant files, and configuration TOMLs for PCGR; tarball and upload to the AWS instance
"""
localrules: pcgr_symlink_somatic, pcgr_symlink_germline, pcgr


import subprocess
from ngs_utils.file_utils import which


rule run_pcgr:
    input:
        vcf = rules.somatic_vcf_filter_pass.output.vcf,
        cns = '{batch}/purple/{batch}.purple.cnv',
        pcgr_data = pcgr_data
    output:
        '{batch}/pcgr/{batch}-somatic.pcgr.html',
        '{batch}/pcgr/{batch}-somatic.pcgr.pass.vcf.gz',
        '{batch}/pcgr/{batch}-somatic.pcgr.snvs_indels.tiers.tsv',
    params:
        output_dir = '{batch}/pcgr',
        genome = run.genome_build,
        sample_name = '{batch}-somatic',
        opt = '--no-docker' if not which('docker') else ''
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
        # TODO: memory based on the mutation number. E.g. over 455k tumor mutations need over 10G
    shell:
        conda_cmd.format('pcgr') +
        'pcgr {input.vcf} {input.cns} -g {params.genome} -o {params.output_dir} -s {params.sample_name} '
        '{params.opt} --pcgr-data {input.pcgr_data}'

rule run_cpsr:
    input:
        vcf = rules.germline_vcf_prep.output.vcf,
        pcgr_data = pcgr_data
    output:
        '{batch}/pcgr/{batch}-normal.cpsr.html'
    params:
        output_dir = '{batch}/pcgr',
        genome_build = run.genome_build,
        sample_name = '{batch}-normal',
        opt = '--no-docker' if not which('docker') else ''
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    shell:
        conda_cmd.format('pcgr') +
        'pcgr {input.vcf} -g {params.genome_build} -o {params.output_dir} -s {params.sample_name} --germline '
        '{params.opt} --pcgr-data {input.pcgr_data} ' \
        '|| echo "Failed to run CPSR" >> {output}'

rule pcgr_symlink_somatic:
    input:
        rules.run_pcgr.output[0]
    output:
        '{batch}/{batch}-somatic.pcgr.html'
    params:
        source = 'pcgr/{batch}-somatic.pcgr.html'
    shell:
        'ln -s {params.source} {output}'

rule pcgr_symlink_germline:
    input:
        rules.run_cpsr.output[0]
    output:
        '{batch}/{batch}-normal.cpsr.html'
    params:
        source = 'pcgr/{batch}-normal.cpsr.html'
    shell:
        'ln -s {params.source} {output}'

######################
###  Target rule

rule pcgr:
    input:
        expand(rules.pcgr_symlink_somatic.output, batch=batch_by_name.keys()),
        expand(rules.pcgr_symlink_germline.output, batch=batch_by_name.keys())
    output:
        temp(touch('log/pcgr.done'))



