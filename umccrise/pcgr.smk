"""
PCGR
-------------
Prepare somatic, germline variant files, and configuration TOMLs for PCGR; tarball and upload to the AWS instance
"""
localrules: pcgr_symlink_somatic, pcgr_symlink_germline, pcgr


from ngs_utils.file_utils import which
from ngs_utils.reference_data import get_predispose_genes_txt, get_predispose_genes_bed

rule run_pcgr:
    input:
        vcf = '{batch}/small_variants/{batch}-somatic.PASS.vcf.gz',
        # cns = '{batch}/purple/{batch}.purple.cnv',
        pcgr_data = hpc.get_ref_file(key='pcgr_data'),
        purple_file = rules.purple_run.output.purity,
    output:
        '{batch}/pcgr/{batch}-somatic.pcgr.html',
        '{batch}/pcgr/{batch}-somatic.pcgr.pass.vcf.gz',
        # '{batch}/pcgr/{batch}-somatic.pcgr.snvs_indels.tiers.tsv',
    params:
        output_dir = '{batch}/pcgr',
        genome = run.genome_build,
        sample_name = '{batch}-somatic',
        opt = '--no-docker' if not which('docker') else ''
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000
        # TODO: memory based on the mutation number. E.g. over 455k tumor mutations need over 10G
    run:
        purity = get_purity(input.purple_file)
        ploidy = get_ploidy(input.purple_file)
        shell(
            conda_cmd.format('pcgr') +
            'pcgr {input.vcf} -g {params.genome} -o {params.output_dir} -s {params.sample_name} '
            '{params.opt} --pcgr-data {input.pcgr_data} --puriry {purity} --ploidy {ploidy}')

rule pcgr_symlink_somatic:
    input:
        rules.run_pcgr.output[0]
    output:
        '{batch}/{batch}-somatic.pcgr.html'
    params:
        source = 'pcgr/{batch}-somatic.pcgr.html'
    shell:
        'ln -s {params.source} {output}'


include_germline = all(b.germline_vcf for b in batch_by_name.values())

if include_germline:
    rule run_cpsr:
        input:
            vcf = 'work/{batch}/small_variants/{batch}-germline.predispose_genes.vcf.gz',
            pcgr_data = hpc.get_ref_file(key='pcgr_data'),
            predispose_bed = get_predispose_genes_bed(run.genome_build),
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
            '{params.opt} --pcgr-data {input.pcgr_data} --predispose-bed {input.predispose_bed}'

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
        expand(rules.pcgr_symlink_germline.output, batch=batch_by_name.keys()) if include_germline else []
    output:
        temp(touch('log/pcgr.done'))

rule cpsr:
    input:
        expand(rules.pcgr_symlink_germline.output, batch=batch_by_name.keys()) if include_germline else []
    output:
        temp(touch('log/cpsr.done'))



