localrules: purple


import glob
import shutil
import platform
from os.path import basename, join
from umccrise import package_path
from reference_data import api as refdata


# localrules: purple
localrules: purple


circos_macos_patch = ('export PERL5LIB=' +
    env_path + '_hmf/lib/site_perl/5.26.2/darwin-thread-multi-2level:' +
    env_path + '_hmf/lib/perl5/site_perl/5.22.0; ') \
    if platform.system() == 'Darwin' \
    else ''



# https://github.com/PapenfussLab/gridss#how-much-memory-should-i-give-gridss
# > At least 4GB + 2GB per thread. It is recommended to run GRIDSS with max heap memory (-Xmx)
# > of 8GB for single-threaded operation (WORKER_THREADS=1), 16GB for multi-core desktop operation,
# > and 31GB for heavily multi-threaded server operation. Note that due to Java's use of Compressed
# > Oops, specifying a max heap size of between 32-48GB effectively reduces the memory available
# > to GRIDSS so is strongly discouraged.
# This should apply to other java apps too
purple_mem = max(min(31000, 4000+2000*threads_per_batch), 8000)
amber_mem  = max(min(31000, 4000+3000*threads_per_batch), 8000)


rule purple_amber:
    input:
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
        het_snps = refdata.get_ref_file(run.genome_build, 'purple_het'),
        ref_fa = refdata.get_ref_file(run.genome_build, 'fa'),
    output:
        'work/{batch}/purple/amber/{batch}.amber.baf.tsv',
        'work/{batch}/purple/amber/{batch}.amber.baf.pcf',
        'work/{batch}/purple/amber/{batch}.amber.baf.vcf.gz',
        'work/{batch}/purple/amber/{batch}.amber.contamination.tsv',
        'work/{batch}/purple/amber/{batch}.amber.contamination.vcf.gz',
        'work/{batch}/purple/amber/{batch}.amber.qc',
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
        normal_name = lambda wc: batch_by_name[wc.batch].normal.name,
        outdir = 'work/{batch}/purple/amber',
        xms = 4000,
    log:
        'log/purple/{batch}/{batch}.amber.log',
    benchmark:
        'benchmarks/{batch}/purple/{batch}-amber.tsv'
    resources:
        mem_mb = int(amber_mem * 1.1),
    threads:
        threads_per_batch
    run:
        shell(
            conda_cmd.format('hmf') +
            f'AMBER -Xms{params.xms}m -Xmx{amber_mem}m '
            f'-tumor {wildcards.batch} '
            f'-tumor_bam {input.tumor_bam} '
            f'-reference {params.normal_name} '
            f'-reference_bam {input.normal_bam} '
            f'-ref_genome {input.ref_fa} '
            f'-loci {input.het_snps} '
            f'-threads {threads} '
            f'-output_dir {params.outdir} 2>&1 | tee {log} '
        )

rule purple_cobalt:
    input:
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
        gc = refdata.get_ref_file(run.genome_build, 'purple_gc'),
        ref_fa = refdata.get_ref_file(run.genome_build, 'fa'),
    output:
        'work/{batch}/purple/cobalt/{batch}.chr.len',
        'work/{batch}/purple/cobalt/{batch}.cobalt.ratio.tsv',
        'work/{batch}/purple/cobalt/{batch}.cobalt.gc.median',
    params:
        outdir = 'work/{batch}/purple/cobalt',
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
        xms = 2000,
    log:
        'log/purple/{batch}/{batch}.cobalt.log'
    benchmark:
        'benchmarks/{batch}/purple/{batch}-cobalt.tsv'
    threads:
        threads_per_batch
    resources:
        mem_mb = int(purple_mem * 1.1),
    run:
        shell(
            conda_cmd.format('hmf') +
            f'COBALT -Xms{params.xms}m -Xmx{purple_mem}m '
            f'-reference {params.normal_sname} '
            f'-reference_bam {input.normal_bam} '
            f'-tumor {wildcards.batch} '
            f'-tumor_bam {input.tumor_bam} '
            f'-ref_genome {input.ref_fa} '
            f'-threads {threads} '
            f'-gc_profile {input.gc} '
            f'-output_dir {params.outdir} 2>&1 | tee {log} '
        )

rule purple_somatic_vcf:
    input:
        '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz',
    output:
        'work/{batch}/purple/somatic.vcf',
    params:
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.rgid,
    group: 'purple_main'
    shell:
        'bcftools view -s {params.tumor_sname} {input} | '
        'bcftools reheader --samples <(echo {wildcards.batch}) > {output}'

rule purple_run:
    input:
        amber_dummy       = 'work/{batch}/purple/amber/{batch}.amber.baf.tsv',
        amber_dummy_pcf   = 'work/{batch}/purple/amber/{batch}.amber.baf.pcf',
        amber_dummy_vcf   = 'work/{batch}/purple/amber/{batch}.amber.baf.vcf.gz',
        cobalt_dummy      = 'work/{batch}/purple/cobalt/{batch}.chr.len',
        cobalt_dummy_pcf  = 'work/{batch}/purple/cobalt/{batch}.cobalt.ratio.tsv',
        manta_sv_filtered = lambda wc: rules.filter_sv_vcf.output.vcf if batch_by_name[wc.batch].sv_vcf else [],
        gc                = refdata.get_ref_file(run.genome_build, 'purple_gc'),
        somatic_vcf       = rules.purple_somatic_vcf.output,
        ref_fa            = refdata.get_ref_file(run.genome_build, 'fa'),
    output:
        cnv           = 'work/{batch}/purple/{batch}.purple.cnv.somatic.tsv',
        gene_cnv      = 'work/{batch}/purple/{batch}.purple.cnv.gene.tsv',
        purity        = 'work/{batch}/purple/{batch}.purple.purity.tsv',
        qc            = 'work/{batch}/purple/{batch}.purple.qc',

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
    log:
        'log/purple/{batch}/{batch}.purple.log'
    benchmark:
        'benchmarks/{batch}/purple/{batch}-purple.tsv'
    threads:
        min(4, threads_per_batch)
    resources:
        mem_mb = int(purple_mem * 1.1),
    run:
        shell(
            conda_cmd.format('hmf') + \
            circos_macos_patch + \
            f'circos -modules ; circos -v ; '
            f'PURPLE -Xms{params.xms}m -Xmx{purple_mem}m '
            f'-amber {params.outdir}/amber '
            f'-cobalt {params.outdir}/cobalt '
            f'-output_dir {params.outdir} '
            f'-reference {params.normal_sname} '
            f'-tumor {wildcards.batch} '
            f'-threads {threads} '
            f'-gc_profile {input.gc} '
            f'{f"-structural_vcf {input.manta_sv_filtered}" if input.manta_sv_filtered else ""} '
            f'-somatic_vcf {input.somatic_vcf} '
            f'-ref_genome {input.ref_fa} '
            f'-circos circos 2>&1 | tee {log} '
        )

rule purple_circos_baf:
    input:
        baf  = 'work/{batch}/purple/circos/{batch}.baf.circos',
        cnv  = 'work/{batch}/purple/circos/{batch}.cnv.circos',
        map  = 'work/{batch}/purple/circos/{batch}.map.circos',
        link = 'work/{batch}/purple/circos/{batch}.link.circos',
        circos_baf_conf = package_path() + '/rmd_files/misc/circos/circos_baf.conf',
    output:
        png = 'work/{batch}/purple/circos_baf/{batch}.circos_baf.png'
    group: 'purple_main'
    params:
        out_dir = 'work/{batch}/purple/circos_baf',
        gaps_txt_prefix = package_path() + '/rmd_files/misc/circos/gaps',
        genome_build = 'hg19' if run.genome_build in ['GRCh37', 'hg19'] else 'hg38',
    run:
        shell('mkdir -p {params.out_dir}')
        shell('cp {params.gaps_txt_prefix}_{params.genome_build}.txt {params.out_dir}/gaps.txt')
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
    group: 'purple_main'
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


#####

rule purple:
    input:
        expand(rules.purple_symlink.output, batch=batch_by_name.keys()),
        expand(rules.purple_circos_baf.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/purple.done'))

