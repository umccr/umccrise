## Cancer gene coverage

localrules: coverage, cacao_symlink_somatic, cacao_symlink_somatic


from ngs_utils.reference_data import get_key_genes_bed
from ngs_utils.file_utils import safe_symlink


min_covs = {'tumor': 30, 'normal': 10}    # minimal coverage to confidently call a variant
max_covs = {'tumor': 200, 'normal': 100}  # coverage above this is suspiciously high (lcr? repeat? CN?)


# Looking at coverage for a limited set of (cancer) genes to assess overall reliability.
# Minimum coverage for normal is 10, 30 for cancer.
#TODO: replace with mosdepth
rule goleft_depth:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        bed = get_key_genes_bed(run.genome_build, coding_only=True),
        ref_fa = ref_fa
    params:
        prefix = lambda wc, output: output[0].replace('.depth.bed', ''),
        cutoff = lambda wc: min_covs[wc.phenotype]
    output:
        '{batch}/coverage/{batch}-{phenotype}.depth.bed'
    threads: threads_per_sample
    resources:
        mem_mb=2000
    shell:
        'goleft depth {input.bam} --reference {ref_fa} --processes {threads} --bed {input.bed} --stats --mincov {params.cutoff} --prefix {params.prefix}'


# Also bringing in global coverage plots for review (tumor only, quick check for CNVs):
rule goleft_plots:
    input:
        bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
    params:
        directory = '{batch}/coverage/{batch}-indexcov',
        xchr = 'X' if run.genome_build == 'GRCh37' else 'chrX'
    output:
        '{batch}/coverage/{batch}-indexcov/index.html'
    resources:
        mem_mb=2000
    shell:
        'goleft indexcov --directory {params.directory} {input.bam} --sex {params.xchr}'


pcgr_genome = 'grch38' if '38' in run.genome_build else 'grch37'

cacao_param = f'--callability_levels_somatic 0:{min_covs["tumor"]}:{max_covs["tumor"]} ' \
              f'--callability_levels_germline 0:{min_covs["normal"]}:{max_covs["normal"]}'

rule run_cacao_somatic:
    input:
        bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
    output:
        report = '{batch}/coverage/cacao_somatic/{batch}.cacao.' + pcgr_genome + '.html'
    params:
        cacao_data = join(loc.extras, 'cacao', 'data'),
        output_dir = '{batch}/coverage/cacao_somatic',
        docker_opt = '--no-docker' if not which('docker') else '',
        output_prefix = '{batch}/coverage/cacao_somatic/{batch}',
    resources:
        mem_mb=2000
    threads: threads_per_sample
    shell:
        conda_cmd.format('pcgr') +
        'cacao_wflow.py {input.bam} {params.cacao_data} {params.output_dir} {pcgr_genome}' \
        ' somatic {params.output_prefix} {params.docker_opt} --threads {threads} ' + cacao_param

rule run_cacao_normal:
    input:
        bam = lambda wc: batch_by_name[wc.batch].normal.bam,
    output:
        report = '{batch}/coverage/cacao_normal/{batch}.cacao.' + pcgr_genome + '.html'
    params:
        cacao_data = join(loc.extras, 'cacao', 'data'),
        output_dir = '{batch}/coverage/cacao_normal',
        docker_opt = '--no-docker' if not which('docker') else '',
        output_prefix = '{batch}/coverage/cacao_somatic/{batch}',
    resources:
        mem_mb=2000
    threads: threads_per_sample
    shell:
        conda_cmd.format('pcgr') +
        'cacao_wflow.py {input.bam} {params.cacao_data} {params.output_dir} {pcgr_genome}'
        ' somatic {params.output_prefix} {params.docker_opt} --threads {threads} ' + cacao_param

rule cacao_symlink_somatic:
    input:
        rules.run_cacao_somatic.output.report
    output:
        '{batch}/{batch}-somatic.cacao.html'
    run:
        safe_symlink(input[0], output[0], rel=True)

rule cacao_symlink_normal:
    input:
        rules.run_cacao_normal.output.report
    output:
        '{batch}/{batch}-normal.cacao.html'
    run:
        safe_symlink(input[0], output[0], rel=True)

rule coverage:
    input:
        expand(rules.goleft_depth.output[0], phenotype=['tumor', 'normal'], batch=batch_by_name.keys()),
        expand(rules.goleft_plots.output[0], batch=batch_by_name.keys()),
        expand(rules.cacao_symlink_somatic.output[0], batch=batch_by_name.keys()),
        expand(rules.cacao_symlink_normal.output[0], batch=batch_by_name.keys()),
    output:
        temp(touch('log/coverage.done'))
