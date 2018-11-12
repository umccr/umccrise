localrules: purple


import glob
import shutil
import platform


rule purple_pileup:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        snp_bed = get_ref_file(run.genome_build, 'purple_het'),
        fasta = ref_fa,
    output:
        'work/{batch}/purple/pileup/{batch}-{phenotype}.mpileup'
    log:
        'log/purple/{batch}/{batch}-{phenotype}.mpileup.log'
    benchmark:
        'benchmarks/{batch}/purple/{batch}-{phenotype}.amber-mpileup.tsv'
    threads:
        max(1, threads_max // len(batch_by_name))
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    shell:
        conda_cmd.format('purple') +
        'sambamba mpileup '
        '-o {output} '
        '-t{threads} '
        '-L <(gunzip -c {input.snp_bed}) '
        '{input.bam} '
        '--samtools -q1 '
        '-f {input.fasta} '
        '2> >(tee -a {log} >&2)'

rule purple_amber:
    input:
        normal_mpileup = 'work/{batch}/purple/pileup/{batch}-normal.mpileup',
        tumor_mpileup  = 'work/{batch}/purple/pileup/{batch}-tumor.mpileup',
    output:
        'work/{batch}/purple/amber/{batch}.amber.baf',
    params:
        outdir = 'work/{batch}/purple/amber',
        jar = join(package_path(), 'amber.jar'),
    log:
        'log/purple/{batch}/{batch}.amber.log',
    shell:
        conda_cmd.format('purple') +
        'java -jar {params.jar} '
        '-sample {wildcards.batch} '
        '-reference {input.normal_mpileup} '
        '-tumor {input.tumor_mpileup} '
        '-output_dir {params.outdir} 2>&1 | tee {log} '

rule purple_cobalt:
    input:
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
        gc = get_ref_file(run.genome_build, 'purple_gc'),
    output:
        'work/{batch}/purple/cobalt/{batch}.cobalt',
    params:
        outdir = 'work/{batch}/purple/cobalt',
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
    log:
        'log/purple/{batch}/{batch}.cobalt.log'
    threads:
        max(1, threads_max // len(batch_by_name))
    shell:
        conda_cmd.format('purple') +
        'COBALT '
        '-reference {params.normal_sname} '
        '-reference_bam {input.normal_bam} '
        '-tumor {wildcards.batch} '
        '-tumor_bam {input.tumor_bam} '
        '-threads {threads} '
        '-gc_profile {input.gc} '
        '-output_dir {params.outdir} 2>&1 | tee {log} '

rule purple_somatic_vcf:
    input:
        rules.somatic_vcf_pon_pass.output.vcf,
    output:
        'work/{batch}/purple/somatic.vcf',
    params:
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
    group: 'purple_run'
    shell:
        'bcftools view -s {params.tumor_sname} {input} | '
        'bcftools reheader --samples <(echo {wildcards.batch}) > {output}'

rule purple_run:
    input:
        cobalt_dummy = 'work/{batch}/purple/cobalt/{batch}.cobalt',
        amber_dummy  = 'work/{batch}/purple/amber/{batch}.amber.baf',
        manta_sv_filtered = rules.filter_sv_vcf.output.vcf,
        gc = get_ref_file(run.genome_build, 'purple_gc'),
        somatic_vcf = rules.purple_somatic_vcf.output,
    output:
        'work/{batch}/purple/{batch}.purple.cnv',
        'work/{batch}/purple/plot/{batch}.circos.png',
    params:
        rundir = 'work/{batch}/purple',
        outdir = 'work/{batch}/purple',
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
        macos_patch = 'export PERL5LIB=' \
            '$CONDA_PREFIX/lib/site_perl/5.26.2/darwin-thread-multi-2level:' \
            '$CONDA_PREFIX/lib/perl5/site_perl/5.22.0 && '\
            if platform.system() == 'Darwin' else ''
    group: 'purple_run'
    threads:
        max(1, threads_max // len(batch_by_name))
    log:
        'log/purple/{batch}/{batch}.purple.log'
    conda:
        'envs/purple.yml'
    shell:
        conda_cmd.format('purple') +
        '{params.macos_patch} '
        'PURPLE '
        '-run_dir {params.rundir} '
        '-output_dir {params.outdir} '
        '-ref_sample {params.normal_sname} '
        '-tumor_sample {wildcards.batch} '
        '-threads {threads} '
        '-gc_profile {input.gc} '
        '-structural_vcf {input.manta_sv_filtered} '
        '-somatic_vcf {input.somatic_vcf} '
        '-circos circos 2>&1 | tee {log} '

rule purple_symlink:
    input:
        'work/{batch}/purple/{batch}.purple.cnv',
        'work/{batch}/purple/plot/{batch}.circos.png',
    output:
        '{batch}/purple/{batch}.purple.cnv',
        '{batch}/purple/{batch}.purple.circos.png',
    params:
        tumor_sname = lambda wc: wc.batch,
        purple_outdir = 'work/{batch}/purple',
    group: 'purple_run'
    run:
        for img_fpath in glob.glob(f'{params.purple_outdir}/plot/*.png'):
            new_name = basename(img_fpath).replace(f'{params.tumor_sname}', f'{wildcards.batch}.purple')
            shutil.copy(img_fpath, join(f'{wildcards.batch}/purple', new_name))

        for fpath in glob.glob(f'{params.purple_outdir}/*.purple.*'):
            new_name = basename(fpath).replace(f'{params.tumor_sname}', f'{wildcards.batch}')
            shutil.copy(fpath, join(f'{wildcards.batch}/purple', new_name))

rule purple:
    input:
        expand(rules.purple_symlink.output, batch=batch_by_name.keys())
    output:
        temp(touch('log/purple.done'))
