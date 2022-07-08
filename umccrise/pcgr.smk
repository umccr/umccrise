"""
PCGR
-------------
Prepare somatic, germline variant files, and configuration TOMLs for PCGR; tarball and upload to the AWS instance
"""
localrules: pcgr, cpsr, cpsr_batch


from os.path import dirname
from ngs_utils.file_utils import which
from ngs_utils.reference_data import get_predispose_genes_txt
from reference_data import api as refdata
from umccrise import get_purity, get_ploidy


# Run PCGR given final prioritised SNVs. Use the PURPLE-inferred purity & ploidy if run.
# In the case of hypermutated samples, the first PCGR run
# (from vcf_stuff/filtering/annotate_somatic_vcf.smk) makes
# sure the intergenic variants are skipped via VEP, which also reduces
# the workload for this run.
rule run_pcgr:
    input:
        vcf = '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz',
        # cns = '{batch}/purple/{batch}.purple.cnv',
        pcgr_data = refdata.get_ref_file(genome=run.genome_build, key='pcgr_data'),
        purple_file = 'work/{batch}/purple/{batch}.purple.purity.tsv' if 'purple' in stages else [],
    output:
        html = 'work/{batch}/pcgr/{batch}-somatic.pcgr.html',
        vcf = 'work/{batch}/pcgr/{batch}-somatic.pcgr.pass.vcf.gz',
        tsv = 'work/{batch}/pcgr/{batch}-somatic.pcgr.snvs_indels.tiers.tsv',
    params:
        genome = run.genome_build,
        sample_name = '{batch}-somatic',
        opt = '--no-docker' if not which('docker') else ''
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000
        # TODO: memory based on the mutation number. E.g. over 455k tumor mutations need over 10G
    group: 'pcgr'
    run:
        output_dir = dirname(output[0])
        cmd = (conda_cmd.format('pcgr') +
            'pcgr {input.vcf} -g {params.genome} -o {output_dir} -s {params.sample_name} '
            '{params.opt} --pcgr-data {input.pcgr_data}')
        if input.purple_file:
            purity = get_purity(input.purple_file)
            ploidy = get_ploidy(input.purple_file)
            cmd += ' --purity {purity} --ploidy {ploidy}'
        shell(cmd)

# copy PCGR html and tsv outputs to <um>/<batch>/ directory
rule pcgr_copy_report:
    input:
        html = rules.run_pcgr.output.html,
        tsv = rules.run_pcgr.output.tsv
    output:
        html = '{batch}/{batch}-somatic.pcgr.html',
        tsv = '{batch}/small_variants/{batch}-somatic.pcgr.snvs_indels.tiers.tsv',
    group: 'pcgr'
    shell:
        'cp {input.html} {output.html}; cp {input.tsv} {output.tsv}'

# run CPSR using merged germline/somatic leakage predisposed VCF and
# predisposed genes txt file
rule run_cpsr:
    input:
        vcf = rules.germline_merge_with_leakage.output.vcf \
              if 'somatic' in stages else \
              rules.germline_predispose_subset.output.vcf,
        pcgr_data = refdata.get_ref_file(genome=run.genome_build, key='pcgr_data'),
        predispose_genes_txt = get_predispose_genes_txt(),
    output:
        'work/{batch}/cpsr/{batch}-normal.cpsr.html'
    params:
        genome_build = run.genome_build,
        sample_name = '{batch}-normal',
        opt = '--no-docker' if not which('docker') else ''
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    group: 'cpsr'
    run:
        output_dir = dirname(output[0])
        shell(conda_cmd.format('pcgr') +
            'pcgr {input.vcf} -g {params.genome_build} -o {output_dir} -s {params.sample_name} '
            '--germline {params.opt} --pcgr-data {input.pcgr_data} '
            '--predispose-genes {input.predispose_genes_txt}')

# copy CPSR html report to <um>/<batch>/
rule cpsr_copy_report:
    input:
        rules.run_cpsr.output[0]
    output:
        '{batch}/{batch}-normal.cpsr.html'
    group: 'cpsr'
    shell:
        'cp {input} {output}'


############

rule pcgr:
    input:
        expand(rules.pcgr_copy_report.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/pcgr.done'))

rule cpsr_batch:
    input:
        lambda wc: rules.cpsr_copy_report.output if batch_by_name[wc.batch].germline_vcf or 'haplotype_caller' in stages else []
    output:
        temp(touch('log/cpsr_{batch}.done'))

rule cpsr:
    input:
        expand(rules.cpsr_batch.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/cpsr.done'))

