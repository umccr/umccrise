## Cancer gene coverage

localrules: coverage, cacao_symlink_somatic, cacao_symlink_normal


from ngs_utils.reference_data import get_key_genes_bed
from ngs_utils.file_utils import safe_symlink
from ngs_utils.file_utils import which


min_covs = {'tumor': 30, 'normal': 10}    # minimal coverage to confidently call a variant
max_covs = {'tumor': 200, 'normal': 100}  # coverage above this is suspiciously high (lcr? repeat? CN?)


# Looking at coverage for a limited set of (cancer) genes to assess overall reliability.
# Minimum coverage for normal is 10, 30 for cancer.
rule mosdepth:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        bed = get_key_genes_bed(run.genome_build, coding_only=True),
    output:
        '{batch}/coverage/{batch}-{phenotype}.quantized.bed.gz',
        '{batch}/coverage/{batch}-{phenotype}.regions.bed.gz',
    params:
        prefix = '{batch}/coverage/{batch}-{phenotype}',
        cutoffs = lambda wc: f'0:1:{min_covs[wc.phenotype]}:{max_covs[wc.phenotype]}:',
    shell:
        'export MOSDEPTH_Q0=NO_COVERAGE && '
        'export MOSDEPTH_Q1=LOW_COVERAGE && '
        'export MOSDEPTH_Q2=CALLABLE && '
        'export MOSDEPTH_Q3=HIGH_COVERAGE && '
        'mosdepth {params.prefix} -q {params.cutoffs} --no-per-base {input.bam} --by {input.bed}'


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
        report = '{batch}/coverage/cacao_somatic/{batch}_' + pcgr_genome + '_coverage_cacao.html'
    params:
        cacao_data = cacao_data,
        output_dir = '{batch}/coverage/cacao_somatic',
        docker_opt = '--no-docker' if not which('docker') else '',
        sample_id = '{batch}',
    resources:
        mem_mb=2000
    threads: threads_per_sample
    shell:
        conda_cmd.format('pcgr') +
        'cacao_wflow.py {input.bam} {params.cacao_data} {params.output_dir} {pcgr_genome}' +
        ' somatic {params.sample_id} {params.docker_opt} --threads {threads} ' + cacao_param

rule run_cacao_normal:
    input:
        bam = lambda wc: batch_by_name[wc.batch].normal.bam,
    output:
        report = '{batch}/coverage/cacao_normal/{batch}_' + pcgr_genome + '_coverage_cacao.html'
    params:
        cacao_data = cacao_data,
        output_dir = '{batch}/coverage/cacao_normal',
        docker_opt = '--no-docker' if not which('docker') else '',
        sample_id = '{batch}',
    resources:
        mem_mb=2000
    threads: threads_per_sample
    shell:
        conda_cmd.format('pcgr') +
        'cacao_wflow.py {input.bam} {params.cacao_data} {params.output_dir} {pcgr_genome}' +
        ' somatic {params.sample_id} {params.docker_opt} --threads {threads} ' + cacao_param

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

rule cacao:
    input:
        expand(rules.cacao_symlink_somatic.output[0], batch=batch_by_name.keys()),
        expand(rules.cacao_symlink_normal.output[0], batch=batch_by_name.keys()),

rule coverage:
    input:
        (expand(rules.mosdepth.output[0], phenotype=['tumor', 'normal'], batch=batch_by_name.keys()) if which('mosdepth') else []),
        expand(rules.goleft_plots.output[0], batch=batch_by_name.keys()),
        expand(rules.cacao_symlink_somatic.output[0], batch=batch_by_name.keys()),
        expand(rules.cacao_symlink_normal.output[0], batch=batch_by_name.keys()),
    output:
        temp(touch('log/coverage.done'))
