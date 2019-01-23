from umccrise.multiqc.prep_data import make_report_metadata, multiqc_prep_data


# localrules: multiqc, copy_logs
localrules: prep_multiqc_data, multiqc, copy_logs


rule prep_multiqc_data:
    input:
        bcbio_mq_filelist = join(run.date_dir, 'multiqc/list_files_final.txt'),
        bcbio_mq_yaml     = join(run.date_dir, 'multiqc/multiqc_config.yaml'),
        bcbio_final_dir   = run.final_dir,
        versions          = 'log/' + run.project_name + '-data_versions.csv',
        programs          = 'log/' + run.project_name + '-programs.txt',
        conpair_contam_t  = lambda wc: f'{wc.batch}/conpair/contamination/{batch_by_name[wc.batch].tumor.name}.txt',
        conpair_contam_n  = lambda wc: f'{wc.batch}/conpair/contamination/{batch_by_name[wc.batch].normal.name}.txt',
        conpair_concord   = lambda wc: f'{wc.batch}/conpair/concordance/{batch_by_name[wc.batch].tumor.name}.txt',
    output:
        filelist            = 'work/{batch}/multiqc_data/filelist.txt',
        generated_conf_yaml = 'work/{batch}/multiqc_data/generated_conf.yaml',
        bcbio_conf_yaml     = 'work/{batch}/multiqc_data/bcbio_conf.yaml',
    params:
        data_dir        = 'work/{batch}/multiqc_data',
        tumor_name      = lambda wc: batch_by_name[wc.batch].tumor.name,
        normal_name     = lambda wc: batch_by_name[wc.batch].normal.name,
    # group: 'multiqc'
    run:
        report_base_path = dirname(abspath(f'{wildcards.batch}/{wildcards.batch}-multiqc_report.html'))
        generated_conf, additional_files = make_report_metadata(
            run,
            tumor_sample=params.tumor_name,
            normal_sample=params.normal_name,
            base_dirpath=report_base_path,
            analysis_dir=run.date_dir,
            program_versions_fpath=input.programs,
            data_versions_fpath=input.versions)
        gold_standard_dir = join(package_path(), 'multiqc', 'gold_standard', 'final.subset.renamed')
        with open(join(gold_standard_dir, 'list_files_final.txt')) as f:
            for l in f:
                l = l.strip()
                additional_files.append(join(gold_standard_dir, l))

        additional_files.extend([input.conpair_contam_t, input.conpair_contam_n, input.conpair_concord])

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


rule batch_multiqc:  # {}
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
        shell(f'LC_ALL=$LC_ALL LANG=$LANG multiqc -f -v -o . -l {list_files} -c {input.bcbio_conf_yaml} -c {input.umccrise_conf_yaml}'
              f' -c {input.generated_conf_yaml} --filename {output.html_file}')


## Additional information
# TODO: link it to MultiQC
rule copy_logs:  # {}
    input:
        versions = join(run.date_dir, 'data_versions.csv'),
        programs = join(run.date_dir, 'programs.txt'),
        conf_dir = run.config_dir
    output:
        versions = 'log/' + run.project_name + '-data_versions.csv',
        programs = 'log/' + run.project_name + '-programs.txt',
        conf_dir = directory(join('log/' + run.project_name + '-config'))
    shell:
        'cp -r {input.versions} {output.versions} && ' \
        'cp -r {input.programs} {output.programs} && ' \
        'cp -r {input.conf_dir} {output.conf_dir}'


rule multiqc:
    input:
        expand(rules.batch_multiqc.output, batch=batch_by_name.keys()),
        rules.copy_logs.output,
    output:
        temp(touch('log/multiqc.done'))












