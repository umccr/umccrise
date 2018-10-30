
rule amber_pileup:
    input:
        bam = lambda wc: bam_from_alias(config, wc.batch, wc.alias) + '.bam'
    output:
        mpileup = join(config['tools']['purple']['outdir'], '{batch}', 'amber/{alias}.mpileup')
    params:
        snp_bed = config['tools']['purple']['hmf_data']['snp_bed'],
        fasta = config['HPC']['ref_fasta']
    log:
        log = join(config['woof']['final_dir'], 'logs', '{batch}/{alias}_amber-pileup.log')
    threads: 32
    shell:
        'echo "[$(date)] start {rule} with wildcards: {wildcards}" > {log.log}; '
        'sambamba mpileup '
        '-t {threads} '
        '-L {params.snp_bed} '
        '{input.bam} '
        '--samtools -q 1 '
        '-f {params.fasta} '
        '> {output.mpileup} 2>> {log.log}; '
        'echo "[$(date)] end {rule} with wildcards: {wildcards}" >> {log.log}; '

rule amber_run:
    input:
        normal_mpileup = lambda wc: join(config['tools']['purple']['outdir'], wc.batch, 'amber', alias_from_pheno(config, wc.batch, 'normal') + '.mpileup'),
        tumor_mpileup = lambda wc: join(config['tools']['purple']['outdir'], wc.batch, 'amber', alias_from_pheno(config, wc.batch, 'tumor') + '.mpileup')
    output:
        tumor_amber = join(config['tools']['purple']['outdir'], '{batch}', 'amber/{tumor_alias}.amber.baf')
    params:
        tumor_alias = lambda wc: alias_from_pheno(config, wc.batch, 'tumor'),
        outdir = join(config['tools']['purple']['outdir'], '{batch}', 'amber'),
        jar = config['tools']['purple']['amber']['jar']
    log:
        log = join(config['woof']['final_dir'], 'logs', '{batch}/{tumor_alias}_amber.log')
    shell:
        'echo "[$(date)] start {rule} with wildcards: {wildcards}" > {log.log}; '
        'java -jar {params.jar} '
        '-sample {params.tumor_alias} '
        '-reference {input.normal_mpileup} '
        '-tumor {input.tumor_mpileup} '
        '-output_dir {params.outdir} >> {log.log} 2>&1; '
        'echo "[$(date)] end {rule} with wildcards: {wildcards}" >> {log.log}; '

rule cobalt_run:
    input:
        normal_bam = lambda wc: bam_from_pheno(config, wc.batch, 'normal') + '.bam',
        tumor_bam = lambda wc: bam_from_pheno(config, wc.batch, 'tumor') + '.bam'
    output:
        tumor_cobalt = join(config['tools']['purple']['outdir'], '{batch}', 'cobalt/{tumor_alias}.cobalt')
    params:
        normal_alias = lambda wc: alias_from_pheno(config, wc.batch, 'normal'),
        tumor_alias = lambda wc: alias_from_pheno(config, wc.batch, 'tumor'),
        gc = config['tools']['purple']['hmf_data']['gc_profile'],
        outdir = join(config['tools']['purple']['outdir'], '{batch}', 'cobalt'),
        jar = config['tools']['purple']['cobalt']['jar']
    log:
        log = join(config['woof']['final_dir'], 'logs', '{batch}/{tumor_alias}_cobalt.log')
    threads: 32
    shell:
        'echo "[$(date)] start {rule} with wildcards: {wildcards}" > {log.log}; '
        'java -jar {params.jar} '
        '-reference {params.normal_alias} '
        '-reference_bam {input.normal_bam} '
        '-tumor {params.tumor_alias} '
        '-tumor_bam {input.tumor_bam} '
        '-threads {threads} '
        '-gc_profile {params.gc} '
        '-output_dir {params.outdir} >> {log.log} 2>&1; '
        'echo "[$(date)] end {rule} with wildcards: {wildcards}" >> {log.log}; '

rule purple_run:
    input:
        cobalt_dummy = lambda wc: join(config['tools']['purple']['outdir'], wc.batch, 'cobalt', alias_from_pheno(config, wc.batch, 'tumor') + '.cobalt'),
        amber_dummy = lambda wc: join(config['tools']['purple']['outdir'], wc.batch, 'amber', alias_from_pheno(config, wc.batch, 'tumor') + '.amber.baf'),
        manta_sv_filtered = lambda wc: join(config['tools']['purple']['outdir'], wc.batch, 'purple', alias_from_pheno(config, wc.batch, 'tumor') + '.manta_filtered.vcf')
    output:
        segs = join(config['tools']['purple']['outdir'], '{batch}', 'purple/{tumor_alias}.purple.cnv')
    params:
        rundir = join(config['tools']['purple']['outdir'], '{batch}'),
        outdir = join(config['tools']['purple']['outdir'], '{batch}', 'purple'),
        jar = config['tools']['purple']['purple']['jar'],
        tumor_alias = lambda wc: alias_from_pheno(config, wc.batch, 'tumor'),
        normal_alias = lambda wc: alias_from_pheno(config, wc.batch, 'normal'),
        gc = config['tools']['purple']['hmf_data']['gc_profile'],
        ensemble_snv = lambda wc: config['bcbio'][wc.batch]['ensemble']
    threads:
        2
    log:
        log = join(config['woof']['final_dir'], 'logs', '{batch}/{tumor_alias}_purple.log')
    shell:
        'echo "[$(date)] start {rule} with wildcards: {wildcards}" > {log.log}; '
        'circos_path=$(which circos); '
        'java -jar {params.jar} '
        '-run_dir {params.rundir} '
        '-output_dir {params.outdir} '
        '-ref_sample {params.normal_alias} '
        '-tumor_sample {params.tumor_alias} '
        '-threads {threads} '
        '-gc_profile {params.gc} '
        '-structural_vcf {input.manta_sv_filtered} '
        '-somatic_vcf {params.ensemble_snv} '
        '-circos ${{circos_path}} >> {log.log} 2>&1; '
        'echo "[$(date)] end {rule} with wildcards: {wildcards}" >> {log.log}; '

