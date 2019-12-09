localrules: purple


import glob
import shutil
import platform
from umccrise import package_path


# localrules: purple
localrules: purple, purple_symlink


circos_macos_patch = ('export PERL5LIB=' +
    env_path + '_hmf/lib/site_perl/5.26.2/darwin-thread-multi-2level:' +
    env_path + '_hmf/lib/perl5/site_perl/5.22.0; ') \
    if platform.system() == 'Darwin' \
    else ''

purple_cpu = min(threads_per_batch, 15)
purple_mem = min(30000, 5000*threads_per_batch)


def get_purple_metric(purple_file, metric='purity'):
    """ Reading the value from somatic sample from Purple output
    """
    with open(purple_file) as f:
        header, values = f.read().split('\n')[:2]
    # #Purity  NormFactor  Score   DiploidProportion  Ploidy  Gender  Status  PolyclonalProportion  MinPurity  MaxPurity  MinPloidy  MaxPloidy  MinDiploidProportion  MaxDiploidProportion  Version  SomaticDeviation
    # 0.7200   1.0400      0.3027  0.8413             1.8611  FEMALE  NORMAL  0.0000                0.6600     0.7700     1.8508     1.8765     0.8241                0.8558                2.17     0.0006
    data = dict(zip(header.strip('#').split('\t'), values.split('\t')))
    purity = float(data[metric])
    return purity


def get_purity(purple_file, phenotype='tumor'):
    """ Reading purity from somatic sample from Purple output
        Assuming purity 100% for normal
    """
    purity = 1.0
    if phenotype == 'tumor':
        purity = get_purple_metric(purple_file, 'purity')
    purity = min(purity, 1.0)
    return purity


def get_ploidy(purple_file):
    return get_purple_metric(purple_file, metric='ploidy')


rule purple_amber:
    input:
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
        snp_bed = hpc.get_ref_file(run.genome_build, 'purple_het'),
        ref_fa = hpc.get_ref_file(run.genome_build, 'fa'),
    output:
        'work/{batch}/purple/amber/{batch}.amber.baf.tsv',
        'work/{batch}/purple/amber/{batch}.amber.baf.vcf.gz',
        'work/{batch}/purple/amber/{batch}.amber.contamination.tsv',
        'work/{batch}/purple/amber/{batch}.amber.contamination.vcf.gz',
        'work/{batch}/purple/amber/{batch}.amber.qc',
        version_build = 'work/{batch}/purple/amber/amber.version',
    params:
        normal_name = lambda wc: batch_by_name[wc.batch].normal.name,
        outdir = 'work/{batch}/purple/amber',
        xms = 5000,
        xmx = purple_mem,
    log:
        'log/purple/{batch}/{batch}.amber.log',
    benchmark:
        'benchmarks/{batch}/purple/{batch}-amber.tsv'
    resources:
        mem_mb = lambda wildcards, attempt: purple_mem + (20000 * attempt),
    threads:
        threads_per_batch
    shell:
        conda_cmd.format('hmf') +
        'AMBER -Xms{params.xms}m -Xmx{params.xmx}m '
        '-tumor {wildcards.batch} '
        '-tumor_bam {input.tumor_bam} '
        '-reference {params.normal_name} '
        '-reference_bam {input.normal_bam} '
        '-ref_genome {input.ref_fa} '
        '-bed {input.snp_bed} '
        '-threads {threads} '
        '-output_dir {params.outdir} 2>&1 | tee {log} '

rule purple_cobalt:
    input:
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
        gc = hpc.get_ref_file(run.genome_build, 'purple_gc'),
        ref_fa = hpc.get_ref_file(run.genome_build, 'fa'),
    output:
        'work/{batch}/purple/cobalt/{batch}.chr.len',
        'work/{batch}/purple/cobalt/{batch}.cobalt.ratio.tsv',
        'work/{batch}/purple/cobalt/{batch}.cobalt.gc.median',
        version_build = 'work/{batch}/purple/cobalt/cobalt.version',
    params:
        outdir = 'work/{batch}/purple/cobalt',
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
        xms = 2000,
        xmx = purple_mem
    log:
        'log/purple/{batch}/{batch}.cobalt.log'
    benchmark:
        'benchmarks/{batch}/purple/{batch}-cobalt.tsv'
    threads:
        threads_per_batch
    resources:
        mem_mb = purple_mem
    shell:
        conda_cmd.format('hmf') +
        'COBALT -Xms{params.xms}m -Xmx{params.xmx}m '
        '-reference {params.normal_sname} '
        '-reference_bam {input.normal_bam} '
        '-tumor {wildcards.batch} '
        '-tumor_bam {input.tumor_bam} '
        '-ref_genome {input.ref_fa} '
        '-threads {threads} '
        '-gc_profile {input.gc} '
        '-output_dir {params.outdir} 2>&1 | tee {log} '

rule purple_somatic_vcf:
    input:
        '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS.vcf.gz',
    output:
        'work/{batch}/purple/somatic.vcf',
    params:
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
    group: 'purple_main'
    shell:
        'bcftools view -s {params.tumor_sname} {input} | '
        'bcftools reheader --samples <(echo {wildcards.batch}) > {output}'

rule purple_run:
    input:
        amber_dummy       = 'work/{batch}/purple/amber/{batch}.amber.baf.tsv',
        amber_dummy_pcf   = 'work/{batch}/purple/amber/{batch}.amber.baf.vcf.gz',
        cobalt_dummy      = 'work/{batch}/purple/cobalt/{batch}.chr.len',
        cobalt_dummy_pcf  = 'work/{batch}/purple/cobalt/{batch}.cobalt.ratio.tsv',
        manta_sv_filtered = rules.filter_sv_vcf.output.vcf,
        gc                = hpc.get_ref_file(run.genome_build, 'purple_gc'),
        somatic_vcf       = rules.purple_somatic_vcf.output,
        ref_fa            = hpc.get_ref_file(run.genome_build, 'fa'),
    output:
        cnv           = 'work/{batch}/purple/{batch}.purple.cnv.somatic.tsv',
        gene_cnv      = 'work/{batch}/purple/{batch}.purple.cnv.gene.tsv',
        purity        = 'work/{batch}/purple/{batch}.purple.purity.tsv',
        qc            = 'work/{batch}/purple/{batch}.purple.qc',
        resc_sv_vcf   = ('work/{batch}/purple/{batch}.purple.sv.vcf.gz' if not is_ffpe else []),
        version_build = 'work/{batch}/purple/purple.version',

        circos_png    = 'work/{batch}/purple/plot/{batch}.circos.png',
        input_png     = 'work/{batch}/purple/plot/{batch}.input.png',
        cn_png        = 'work/{batch}/purple/plot/{batch}.copynumber.png',
        ma_png        = 'work/{batch}/purple/plot/{batch}.map.png',
        purity_png    = 'work/{batch}/purple/plot/{batch}.purity.range.png',
        segment_png   = 'work/{batch}/purple/plot/{batch}.segment.png',
        clonality_png = 'work/{batch}/purple/plot/{batch}.somatic.clonality.png',
        ploidy_png    = 'work/{batch}/purple/plot/{batch}.somatic.png',
        rainfall_png  = 'work/{batch}/purple/plot/{batch}.somatic.rainfall.png',

        baf           = 'work/{batch}/purple/circos/{batch}.baf.circos',
        cnv_circos    = 'work/{batch}/purple/circos/{batch}.cnv.circos',
        map           = 'work/{batch}/purple/circos/{batch}.map.circos',
        link          = 'work/{batch}/purple/circos/{batch}.link.circos',
    group: 'purple_main'
    params:
        outdir = 'work/{batch}/purple',
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
        xms = 2000,
        xmx = purple_mem,
    log:
        'log/purple/{batch}/{batch}.purple.log'
    benchmark:
        'benchmarks/{batch}/purple/{batch}-purple.tsv'
    threads:
        threads_per_batch
    resources:
        mem_mb = purple_mem
    shell:
       conda_cmd.format('hmf') + \
        circos_macos_patch + \
        'circos -modules ; circos -v ; '
        'PURPLE -Xms{params.xms}m -Xmx{params.xmx}m '
        '-amber {params.outdir}/amber '
        '-cobalt {params.outdir}/cobalt '
        '-output_dir {params.outdir} '
        '-reference {params.normal_sname} '
        '-tumor {wildcards.batch} '
        '-threads {threads} '
        '-gc_profile {input.gc} '
         # '-sv_recovery_vcf {input.manta_sv_filtered} '
        '-structural_vcf {input.manta_sv_filtered} '
        '-somatic_vcf {input.somatic_vcf} '
        '-ref_genome {input.ref_fa} '
        '-circos circos 2>&1 | tee {log} '

rule purple_circos_baf:
    input:
        baf  = 'work/{batch}/purple/circos/{batch}.baf.circos',
        cnv  = 'work/{batch}/purple/circos/{batch}.cnv.circos',
        map  = 'work/{batch}/purple/circos/{batch}.map.circos',
        link = 'work/{batch}/purple/circos/{batch}.link.circos',
        circos_baf_conf = package_path() + '/rmd_files/misc/circos/circos_baf.conf',
    output:
        png = 'work/{batch}/purple/circos_baf/{batch}.circos_baf.png'
    params:
        out_dir = 'work/{batch}/purple/circos_baf',
        gaps_txt_prefix = package_path() + '/rmd_files/misc/circos/gaps',
        genome_build = 'hg19' if run.genome_build in ['GRCh37', 'hg19'] else 'hg38',
    run:
        shell('mkdir -p {params.out_dir}')
        shell('cp {params.gaps_txt_prefix}_{params.genome_build}.txt {params.out_dir}/gaps.txt')
        shell('cp {input.ideo_conf} {params.out_dir}')
        out_conf = join(params.out_dir, basename(input.circos_baf_conf))
        shell('sed s/SAMPLE/{wildcards.batch}/ {input.circos_baf_conf} > ' + out_conf)
        shell('cp {input.baf} {params.out_dir}')
        shell('cp {input.cnv} {params.out_dir}')
        shell('cp {input.map} {params.out_dir}')
        shell('cp {input.link} {params.out_dir}')
        out_file = basename(output.png)
        shell(conda_cmd.format('hmf') +
              'circos -modules ; circos -v ; ' +
              circos_macos_patch +
              'circos -nosvg -conf ' + out_conf + ' -outputdir {params.out_dir} -outputfile {out_file}')

rule purple_symlink:
    input:
        rules.purple_run.output.cnv,
        rules.purple_run.output.gene_cnv,
        rules.purple_run.output.circos_png,
        rules.purple_circos_baf.output.png,
    output:
        cnv            = '{batch}/purple/{batch}.purple.cnv.somatic.tsv',
        gene_cnv       = '{batch}/purple/{batch}.purple.cnv.gene.tsv',
        circos_png     = '{batch}/purple/{batch}.purple.circos.png',
        circos_baf_png = '{batch}/purple/{batch}.purple.circos_baf.png',
    params:
        tumor_sname = lambda wc: wc.batch,
        purple_outdir = 'work/{batch}/purple',
    run:
        for img_fpath in glob.glob(f'{params.purple_outdir}/plot/*.png'):
            new_name = basename(img_fpath).replace(f'{params.tumor_sname}', f'{wildcards.batch}.purple')
            shutil.copy(img_fpath, join(f'{wildcards.batch}/purple', new_name))

        for img_fpath in glob.glob(f'{params.purple_outdir}/circos_baf/*.png'):
            new_name = basename(img_fpath).replace(f'{params.tumor_sname}', f'{wildcards.batch}.purple')
            shutil.copy(img_fpath, join(f'{wildcards.batch}/purple', new_name))

        for fpath in glob.glob(f'{params.purple_outdir}/*.purple.*'):
            new_name = basename(fpath).replace(f'{params.tumor_sname}', f'{wildcards.batch}')
            shutil.copy(fpath, join(f'{wildcards.batch}/purple', new_name))


rule purple:
    input:
        expand(rules.purple_symlink.output, batch=batch_by_name.keys()),
        expand(rules.purple_circos_baf.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/purple.done'))

