localrules: purple


import glob
import shutil
import platform


# localrules: purple
localrules: purple, purple_symlink, purple_somatic_vcf, copy_inputs_from_bcbio


circos_macos_patch = ('export PERL5LIB=' +
    env_path + '_purple/lib/site_perl/5.26.2/darwin-thread-multi-2level:' +
    env_path + '_purple/lib/perl5/site_perl/5.22.0 && ') if platform.system() == 'Darwin' else ''


if glob.glob(join(run.work_dir, f'structural/*/purple/amber')):
    rule copy_inputs_from_bcbio:
        input:
            bcbio_amber_dir  = lambda wc: join(run.work_dir, f'structural/{batch_by_name[wc.batch].tumor.name}/purple/amber'),
            bcbio_cobalt_dir = lambda wc: join(run.work_dir, f'structural/{batch_by_name[wc.batch].tumor.name}/purple/cobalt'),
        output:
            'work/{batch}/purple/cobalt/{batch}.cobalt',
            'work/{batch}/purple/cobalt/{batch}.cobalt.ratio.pcf',
            'work/{batch}/purple/amber/{batch}.amber.baf',
        # group: 'purple'
        params:
            tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
            cobalt_dir = 'work/{batch}/purple/cobalt',
            amber_dir = 'work/{batch}/purple/amber',
        run:
            for bcbio_fp in glob.glob(f'{input.bcbio_amber_dir}/*'):
                bcbio_fn = basename(bcbio_fp)
                umccrise_fn = bcbio_fn.replace(params.tumor_name, wildcards.batch)
                umccrise_fp = join(params.amber_dir, umccrise_fn)
                shell(f'cp {bcbio_fp} {umccrise_fp}')
            for bcbio_fp in glob.glob(f'{input.bcbio_cobalt_dir}/*'):
                bcbio_fn = basename(bcbio_fp)
                umccrise_fn = bcbio_fn.replace(params.tumor_name, wildcards.batch)
                umccrise_fp = join(params.cobalt_dir, umccrise_fn)
                shell(f'cp {bcbio_fp} {umccrise_fp}')

else:
    rule purple_amber:
        input:
            tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
            normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
            snp_bed = get_ref_file(run.genome_build, 'purple_het'),
            ref_fa = ref_fa,
        output:
            'work/{batch}/purple/amber/{batch}.amber.baf',
        # group: 'purple'
        params:
            normal_name = lambda wc: batch_by_name[wc.batch].normal.name,
            outdir = 'work/{batch}/purple/amber',
            jar = join(package_path(), 'jars', 'amber.jar'),
            xms = 5000,
            xmx = min(40000, 8000*threads_per_batch),
        log:
            'log/purple/{batch}/{batch}.amber.log',
        benchmark:
            'benchmarks/{batch}/purple/{batch}-amber.tsv'
        resources:
            mem_mb = min(50000, 10000*threads_per_batch),
        threads:
            threads_per_batch
        shell:
            conda_cmd.format('purple') +
            'java -Xms{params.xms}m -Xmx{params.xmx}m -jar {params.jar} '
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
            gc = get_ref_file(run.genome_build, 'purple_gc'),
        output:
            'work/{batch}/purple/cobalt/{batch}.cobalt',
            'work/{batch}/purple/cobalt/{batch}.cobalt.ratio.pcf',
        # group: 'purple'
        params:
            outdir = 'work/{batch}/purple/cobalt',
            normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
            xms = 2000,
            xmx = min(40000, 3000*threads_per_batch),
        log:
            'log/purple/{batch}/{batch}.cobalt.log'
        benchmark:
            'benchmarks/{batch}/purple/{batch}-cobalt.tsv'
        threads:
            threads_per_batch
        resources:
            mem_mb = min(50000, 3500*threads_per_batch)
        shell:
            conda_cmd.format('purple') +
            'COBALT -Xms{params.xms}m -Xmx{params.xmx}m '
            '-reference {params.normal_sname} '
            '-reference_bam {input.normal_bam} '
            '-tumor {wildcards.batch} '
            '-tumor_bam {input.tumor_bam} '
            '-threads {threads} '
            '-gc_profile {input.gc} '
            '-output_dir {params.outdir} 2>&1 | tee {log} '

rule purple_somatic_vcf:
    input:
        rules.somatic_vcf_filter_pass.output.vcf,
    output:
        'work/{batch}/purple/somatic.vcf',
    params:
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
    # group: 'purple'
    shell:
        'bcftools view -s {params.tumor_sname} {input} | '
        'bcftools reheader --samples <(echo {wildcards.batch}) > {output}'

rule purple_run:
    input:
        cobalt_dummy      = 'work/{batch}/purple/cobalt/{batch}.cobalt',
        cobalt_dummy_pcf  = 'work/{batch}/purple/cobalt/{batch}.cobalt.ratio.pcf',
        amber_dummy       = 'work/{batch}/purple/amber/{batch}.amber.baf',
        manta_sv_filtered = rules.filter_sv_vcf.output.vcf,
        gc                = get_ref_file(run.genome_build, 'purple_gc'),
        somatic_vcf       = rules.purple_somatic_vcf.output,
    output:
        cnv          = 'work/{batch}/purple/{batch}.purple.cnv',
        gene_cnv     = 'work/{batch}/purple/{batch}.purple.gene.cnv',
        germline_cnv = 'work/{batch}/purple/{batch}.purple.germline.cnv',
        circos_png   = 'work/{batch}/purple/plot/{batch}.circos.png',
        input_png    = 'work/{batch}/purple/plot/{batch}.input.png',
        cn_png       = 'work/{batch}/purple/plot/{batch}.copyNumber.png',
        ma_png       = 'work/{batch}/purple/plot/{batch}.minor_allele.png',
        variant_png  = 'work/{batch}/purple/plot/{batch}.variant.png',
        purity       = 'work/{batch}/purple/{batch}.purple.purity',
        baf          = 'work/{batch}/purple/circos/{batch}.baf.circos',
        cnv_circos   = 'work/{batch}/purple/circos/{batch}.cnv.circos',
        map          = 'work/{batch}/purple/circos/{batch}.map.circos',
        link         = 'work/{batch}/purple/circos/{batch}.link.circos',
    # group: 'purple'
    params:
        jar = join(package_path(), 'jars', 'purple.jar'),
        rundir = 'work/{batch}/purple',
        outdir = 'work/{batch}/purple',
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
        xms = 2000,
        xmx = min(40000, 3000*threads_per_batch),
    log:
        'log/purple/{batch}/{batch}.purple.log'
    benchmark:
        'benchmarks/{batch}/purple/{batch}-purple.tsv'
    threads:
        threads_per_batch
    resources:
        mem_mb = min(50000, 3500*threads_per_batch)
    shell:
        conda_cmd.format('purple') +
        'circos -modules ; circos -v ; ' +
        circos_macos_patch +
        'java -Xms{params.xms}m -Xmx{params.xmx}m -jar {params.jar} '
        '-run_dir {params.rundir} '
        '-output_dir {params.outdir} '
        '-ref_sample {params.normal_sname} '
        '-tumor_sample {wildcards.batch} '
        '-threads {threads} '
        '-gc_profile {input.gc} '
        '-structural_vcf {input.manta_sv_filtered} '
        '-somatic_vcf {input.somatic_vcf} '
        '-circos circos 2>&1 | tee {log} '

rule purple_circos_baf:
    input:
        baf = 'work/{batch}/purple/circos/{batch}.baf.circos',
        cnv = 'work/{batch}/purple/circos/{batch}.cnv.circos',
        map = 'work/{batch}/purple/circos/{batch}.map.circos',
        link = 'work/{batch}/purple/circos/{batch}.link.circos',
        circos_baf_conf = package_path() + '/rmd_files/templates/circos/circos_baf.conf',
        gaps_txt = package_path() + '/rmd_files/templates/circos/gaps.txt',
        ideo_conf = package_path() + '/rmd_files/templates/circos/ideogram.conf',
    output:
        png = 'work/{batch}/purple/circos_baf/{batch}.circos_baf.png'
    # group: 'purple'
    params:
        out_dir = 'work/{batch}/purple/circos_baf'
    run:
        shell('mkdir -p {params.out_dir}')
        shell('cp {input.gaps_txt} {params.out_dir}')
        shell('cp {input.ideo_conf} {params.out_dir}')
        out_conf = join(params.out_dir, basename(input.circos_baf_conf))
        shell('sed s/SAMPLE/{wildcards.batch}/ {input.circos_baf_conf} > ' + out_conf)
        shell('cp {input.baf} {params.out_dir}')
        shell('cp {input.cnv} {params.out_dir}')
        shell('cp {input.map} {params.out_dir}')
        shell('cp {input.link} {params.out_dir}')
        out_file = basename(output.png)
        shell(conda_cmd.format('purple') +
              'circos -modules ; circos -v ; ' +
              circos_macos_patch +
              'circos -nosvg -conf ' + out_conf + ' -outputdir {params.out_dir} -outputfile {out_file}')

rule purple_symlink:
    input:
        rules.purple_run.output.cnv,
        rules.purple_run.output.gene_cnv,
        rules.purple_run.output.circos_png,
    output:
        cnv = '{batch}/purple/{batch}.purple.cnv',
        gene_cnv = '{batch}/purple/{batch}.purple.gene.cnv',
        circos_png = '{batch}/purple/{batch}.purple.circos.png',
    params:
        tumor_sname = lambda wc: wc.batch,
        purple_outdir = 'work/{batch}/purple',
    # group: 'purple'
    run:
        for img_fpath in glob.glob(f'{params.purple_outdir}/plot/*.png'):
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

