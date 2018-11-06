import glob
from ngs_utils.file_utils import symlink_plus, safe_symlink

rule amber_pileup:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        snp_bed = join(loc.purple, 'GermlineHetPon.' + run.genome_build + '.bed'),
        fasta = ref_fa,
    output:
        mpileup = '{batch}/purple/pileup/{batch}-{phenotype}.mpileup'
    log:
        'log/purple/{batch}/{batch}-{phenotype}.mpileup.log'
    benchmark:
        '{batch}/purple/benchmarks/{batch}-{phenotype}.amber-mpileup.tsv'
    threads: max(1, threads_max // len(batch_by_name))
    shell:
        'sambamba mpileup '
        '-t {threads} '
        '-L {input.snp_bed} '
        '{input.bam} '
        '--samtools -q 1 '
        '-f {input.fasta} '
        '> {output.mpileup} 2>> {log}; '

rule amber_run:
    input:
        normal_mpileup = '{batch}/purple/pileup/{batch}-normal.mpileup',
        tumor_mpileup  = '{batch}/purple/pileup/{batch}-tumor.mpileup',
    output:
        directory('{batch}/purple/amber'),
    params:
        outdir = '{batch}/purple/amber',
        jar = join(loc.purple, 'amber-1.6.jar'),
        tumor_sname = lambda wc: batch_by_name[wc.batch].tumor.name,
    log:
        'log/purple/{batch}/{batch}.amber.log',
    shell:
        'java -jar {params.jar} '
        '-sample {params.tumor_sname} '
        '-reference {input.normal_mpileup} '
        '-tumor {input.tumor_mpileup} '
        '-output_dir {params.outdir} 2>&1 | tee {log} '

rule cobalt_run:
    input:
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
        gc = join(loc.purple, 'GC_profile.' + run.genome_build + '.1000bp.cnp'),
    output:
        directory('{batch}/purple/cobalt'),
    params:
        jar = join(loc.purple, 'cobalt-1.4.jar'),
        outdir = '{batch}/purple/cobalt',
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
    log:
        'log/purple/{batch}/{batch}.cobalt.log'
    threads: max(1, threads_max // len(batch_by_name))
    shell:
        'java -jar {params.jar} '
        '-reference {params.normal_sname} '
        '-reference_bam {input.normal_bam} '
        '-tumor {params.tumor_sname} '
        '-tumor_bam {input.tumor_bam} '
        '-threads {threads} '
        '-gc_profile {input.gc} '
        '-output_dir {params.outdir} 2>&1 | tee {log} '

rule purple_run:
    input:
        cobalt = directory('{batch}/purple/cobalt'),
        amber = directory('{batch}/purple/amber'),
        manta_sv_filtered = rules.filter_sv_vcf.output.vcf,
        somatic_vcf = rules.somatic_vcf_pon_pass.output.vcf,
        gc = join(loc.purple, 'GC_profile.' + run.genome_build + '.1000bp.cnp'),
    output:
        directory('{batch}/purple/purple_perse'),
    params:
        rundir = '{batch}/purple',
        outdir = '{batch}/purple/purple_perse',
        jar = join(loc.purple, 'purple-2.15.jar'),
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
    threads:
        max(1, threads_max // len(batch_by_name))
    log:
        'log/purple/{batch}/{batch}.purple.log'
    shell:
        'circos_path=$(which circos); '
        'java -jar {params.jar} '
        '-run_dir {params.rundir} '
        '-output_dir {params.outdir} '
        '-ref_sample {params.normal_sname} '
        '-tumor_sample {params.tumor_sname} '
        '-threads {threads} '
        '-gc_profile {input.gc} '
        '-structural_vcf {input.manta_sv_filtered} '
        '-somatic_vcf {input.somatic_vcf} '
        '-circos ${{circos_path}} 2>&1 | tee {log} '

rule purple_symlink:
    input:
        purple_perse = rules.purple_run.output,
    output:
        '{batch}/purple/{batch}.purple.cnv',
        '{batch}/purple/{batch}.purple.circos.png',
    params:
        tumor_sname = lambda wc: batch_by_name[wc.batch].tumor.name,
        batch = lambda wc: wc.batch,
    run:
        for img_fpath in glob.glob(f'{input.purple_perse}/plot/*.png'):
            new_name = basename(img_fpath).replace(f'{params.tumor_sname}', f'{params.batch}.purple')
            safe_symlink(img_fpath, join(f'{params.batch}/purple', new_name), rel=True)

        for fpath in glob.glob(f'{input.purple_perse}/*.purple.*'):
            new_name = basename(fpath).replace(f'{params.tumor_sname}', f'{params.batch}')
            safe_symlink(fpath, join(f'{params.batch}/purple', new_name), rel=True)

rule purple:
    input:
        expand(rules.purple_symlink.output, batch=batch_by_name.keys())
    output:
        temp(touch('log/purple.done'))
