from umccrise.multiqc.prep_data import make_report_metadata, multiqc_prep_data
import cyvcf2
import yaml


# localrules: multiqc, copy_logs
localrules: prep_multiqc_data, multiqc, copy_config, copy_logs


versions = join(run.date_dir, 'data_versions.csv')
programs = join(run.date_dir, 'programs.txt')


rule copy_config:
    input:
        conf_dir = run.config_dir
    output:
        conf_dir = directory(join('log/config'))
    shell:
        'cp -r {input.conf_dir} {output.conf_dir}'


rule copy_logs:
    input:
        versions = versions,
        programs = programs,
    output:
        versions = 'log/data_versions.csv',
        programs = 'log/programs.txt',
    shell:
        'cp -r {input.versions} {output.versions} && ' \
        'cp -r {input.programs} {output.programs}'


rule varqc_report:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf
    output:
        'work/{batch}/small_variants/pon_mqc.yaml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name,
    run:
        total_cnt = 0
        pass_cnt = 0
        pon_cnt = 0
        pon_filt = 0
        vcf = cyvcf2.VCF(input.vcf)
        for rec in vcf:
            total_cnt += 1
            if rec.FILTER == 'PASS':
                pass_cnt += 1
            if 'PoN' in rec.FILTER:
                pon_filt += 1
            if rec.INFO.get('PoN', 0) > 0:
                pon_cnt += 1
        pass_pct = pass_cnt / total_cnt
        pon_pct = pon_cnt / total_cnt
        pon_filt_pct = pon_filt / total_cnt

        data = dict(
            # id='umccr_pon',
            # section_name='UMCCR panel of normals',
            # description='Percentage of variants found in UMCCR panel of normals',
            plot_type='generalstats',
            pconfig=[{
                'PoN': dict(
                    min=0,
                    max=100,
                    scale='RdYlGn',
                    suffix='%',
                ),
                'PASS': dict(
                    min=0,
                    max=100,
                    scale='RdYlGn',
                    suffix='%',
                ),
                'PoN filt': dict(
                    min=0,
                    max=100,
                    scale='RdYlGn',
                    suffix='%',
                    hidden=True,
                ),
            }],
            data={
                params.sample: {
                    'PoN': pon_pct,
                    'PASS': pass_cnt,
                    'PoN filt': pon_filt_pct
                }
            }
        )
        with open(output[0], 'w') as out:
            yaml.dump(data, out)


rule prep_multiqc_data:
    input:
        bcbio_mq_filelist = join(run.date_dir, 'multiqc/list_files_final.txt'),
        bcbio_mq_yaml     = join(run.date_dir, 'multiqc/multiqc_config.yaml'),
        bcbio_final_dir   = run.final_dir,
        conpair_concord   = rules.run_conpair.output.concord,
        conpair_contam    = rules.run_conpair.output.contam,
        varqc_report      = rules.varqc_report.output[0],
    output:
        filelist            = 'work/{batch}/multiqc_data/filelist.txt',
        generated_conf_yaml = 'work/{batch}/multiqc_data/generated_conf.yaml',
        bcbio_conf_yaml     = 'work/{batch}/multiqc_data/bcbio_conf.yaml',
    params:
        data_dir        = 'work/{batch}/multiqc_data',
        tumor_name      = lambda wc: batch_by_name[wc.batch].tumor.name,
        normal_name     = lambda wc: batch_by_name[wc.batch].normal.name,
        versions        = 'log/data_versions.csv',
        programs        = 'log/programs.txt',
    # group: 'multiqc'
    run:
        report_base_path = dirname(abspath(f'{wildcards.batch}/{wildcards.batch}-multiqc_report.html'))
        generated_conf, additional_files = make_report_metadata(
            run,
            tumor_sample=params.tumor_name,
            normal_sample=params.normal_name,
            base_dirpath=report_base_path,
            analysis_dir=run.date_dir,
            program_versions_fpath=verify_file(params.programs, silent=True),
            data_versions_fpath=verify_file(params.versions, silent=True)
        )
        gold_standard_dir = join(package_path(), 'multiqc', 'gold_standard', 'final.subset.renamed')
        with open(join(gold_standard_dir, 'list_files_final.txt')) as f:
            for l in f:
                l = l.strip()
                additional_files.append(join(gold_standard_dir, l))

        additional_files.extend([
            join(input.conpair_contam, params.tumor_name + '.txt'),
            join(input.conpair_contam, params.normal_name + '.txt'),
            join(input.conpair_concord, params.tumor_name + '.txt'),
            varqc_report,
        ])

        multiqc_prep_data(
            bcbio_mq_filelist=input.bcbio_mq_filelist,
            bcbio_mq_yaml=input.bcbio_mq_yaml,
            bcbio_final_dir=input.bcbio_final_dir,
            new_mq_data_dir=params.data_dir,
            generated_conf=generated_conf,
            out_filelist_file=output.filelist,
            out_conf_yaml=output.generated_conf_yaml,
            new_bcbio_mq_yaml=output.bcbio_conf_yaml,
            additional_files=additional_files,
            exclude_files=r'.*qsignature.*',
        )


rule batch_multiqc:
    input:
        filelist            = 'work/{batch}/multiqc_data/filelist.txt',
        generated_conf_yaml = 'work/{batch}/multiqc_data/generated_conf.yaml',
        bcbio_conf_yaml     = 'work/{batch}/multiqc_data/bcbio_conf.yaml',
        umccrise_conf_yaml  = join(package_path(), 'multiqc', 'multiqc_config.yaml'),
    output:
        html_file           = '{batch}/{batch}-multiqc_report.html'
    # group: 'multiqc'
    run:
        other_samples=[
            s.name for s in run.samples if s.name not in [batch_by_name[wildcards.batch].tumor.name,
                                                          batch_by_name[wildcards.batch].normal.name]]
        if other_samples:
            greps = ''.join(f' | grep -v {sn}' for sn in other_samples)
            list_files = f'<(cat {input.filelist}{greps})'
        else:
            list_files = input.filelist
        shell(f'LC_ALL=$LC_ALL LANG=$LANG multiqc -f -o . -l {list_files}'
              f' -c {input.bcbio_conf_yaml} -c {input.umccrise_conf_yaml}'
              f' -c {input.generated_conf_yaml} --filename {output.html_file}')


rule multiqc:
    input:
        expand(rules.batch_multiqc.output, batch=batch_by_name.keys()),
        rules.copy_config.output,
        (rules.copy_logs.output if all(isfile(fp) for fp in rules.copy_logs.input) else []),
    output:
        temp(touch('log/multiqc.done'))













