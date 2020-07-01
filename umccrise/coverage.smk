import subprocess
from os.path import abspath, dirname, join
from ngs_utils.reference_data import get_key_genes_bed
from ngs_utils.file_utils import safe_symlink
from ngs_utils.file_utils import which
from umccrise import get_purity
from reference_data import api as refdata


localrules: coverage, cacao, goleft, mosdepth


MIN_VD = 12     # minimal coverage to call a pure heterozygous variant
HIGH_VD = 100   # purity-normalized coverage above this is suspiciously high (lcr? repeat? CN?)


def _get_low_high_covs(phenotype, purple_file):
    # pure min cov  purity  min cov
    # 12            0.4     30
    # 12            0.5     24
    # 12            0.6     20
    # 12            0.7     17
    # 12            0.8     15
    # 12            0.9     13
    # 12            1.0     12
    p = get_purity(purple_file, phenotype) if purple_file else 1.0
    return round(MIN_VD / p), round(HIGH_VD / p)


# Looking at coverage for a limited set of (cancer) genes to assess overall reliability.
rule run_mosdepth:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        bed = get_key_genes_bed(run.genome_build, coding_only=True),
        purple_file = rules.purple_run.output.purity if 'purple' in stages else [],
        ref_fa = refdata.get_ref_file(run.genome_build, 'fa'),
    output:
        '{batch}/coverage/{batch}-{phenotype}.quantized.bed.gz',
        '{batch}/coverage/{batch}-{phenotype}.regions.bed.gz',
    params:
        prefix = '{batch}/coverage/{batch}-{phenotype}',
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000
    threads:
        threads_per_sample
    run:
        low_cov, high_cov = _get_low_high_covs(wildcards.phenotype, input.purple_file)
        cutoffs = f'0:1:{low_cov}:{high_cov}:'
        mosdepth_cmd = (
            'export MOSDEPTH_Q0=NO_COVERAGE && '
            'export MOSDEPTH_Q1=LOW_COVERAGE && '
            'export MOSDEPTH_Q2=CALLABLE && '
            'export MOSDEPTH_Q3=HIGH_COVERAGE && '
            f'mosdepth {params.prefix} {input.bam} '
            f'-q {cutoffs} '
            f'--fasta {input.ref_fa} '
            f'--by {input.bed} '
            f'--no-per-base '
        )
        shell(mosdepth_cmd)


# Also bringing in global coverage plots for review (tumor only, quick check for CNVs):
# For CRAM, indexcov takes the .crai files directly on input.
#  crai resolution is often much less than 100KB (compared to) 16KB for the bam index,
#  but it is sufficient to find large-scale differences in coverage.
rule goleft_plots:
    input:
        bam = lambda wc: batch_by_name[wc.batch].tumor.bam + \
            ('.crai' if batch_by_name[wc.batch].tumor.bam.endswith('.cram') else ''),
        fai = refdata.get_ref_file(run.genome_build, 'fa') + '.fai',
    params:
        directory = '{batch}/coverage/{batch}-indexcov',
        xchr = 'X' if run.genome_build == 'GRCh37' else 'chrX'
    output:
        '{batch}/coverage/{batch}-indexcov/index.html'
    resources:
        mem_mb=2000
    shell:
        'goleft indexcov {input.bam} '
        '--fai {input.fai} '
        '--directory {params.directory} '
        '--sex {params.xchr} '
        '2>&1 | '
        'grep -v "indexcov: excluding chromosome" | '
        'grep -v "no reference stats found" >&2'


pcgr_genome = 'grch38' if '38' in run.genome_build else 'grch37'

rule run_cacao:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        purple_file = rules.purple_run.output.purity if 'purple' in stages else [],
        ref_fa = refdata.get_ref_file(run.genome_build, 'fa'),
    output:
        report = '{batch}/coverage/cacao_{phenotype}/{batch}_' + pcgr_genome + '_coverage_cacao.html'
    params:
        cacao_data = refdata.get_ref_file(genome=run.genome_build, key='cacao_data'),
        output_dir = '{batch}/coverage/cacao_{phenotype}',
        docker_opt = '--no-docker' if not which('docker') else '',
        sample_id = '{batch}',
        mode = lambda wc: 'hereditary' if wc.phenotype == 'normal' else 'somatic',
        levels = lambda wc: 'germline' if wc.phenotype == 'normal' else 'somatic',
    resources:
        mem_mb=8000
    threads: threads_per_sample
    group: 'cacao'
    run:
        low_cov, high_cov = _get_low_high_covs(wildcards.phenotype, input.purple_file)
        cutoffs = f'0:{low_cov}:{high_cov}'
        shell(
            conda_cmd.format('pcgr') +
            f'cacao_wflow.py {input.bam} {params.cacao_data} {params.output_dir} {pcgr_genome}' +
            f' {params.mode} {params.sample_id} {params.docker_opt} --threads {threads}'
            f' --callability_levels_{params.levels} {cutoffs} --ref-fasta {input.ref_fa}'
        )

rule cacao_symlink:
    input:
        rules.run_cacao.output.report
    output:
        '{batch}/{batch}-{phenotype}.cacao.html'
    group: 'cacao'
    run:
        safe_symlink(input[0], output[0], rel=True)


#####

rule cacao:
    input:
        expand(rules.cacao_symlink.output[0], batch=batch_by_name.keys(), phenotype=['tumor', 'normal'])
    output:
        temp(touch('log/coverage.done'))


rule mosdepth:
    input:
        expand(rules.run_mosdepth.output[0], batch=batch_by_name.keys(), phenotype=['tumor', 'normal'])
    output:
        temp(touch('log/mosdepth.done'))


rule goleft:
    input:
        expand(rules.goleft_plots.output[0], batch=batch_by_name.keys())
    output:
        temp(touch('log/goleft.done'))


rule coverage:
    input:
        'log/coverage.done',
        'log/mosdepth.done',
        'log/goleft.done',
    output:
        temp(touch('log/coverage.done'))


