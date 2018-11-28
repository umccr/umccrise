## IGV

localrules: multiqc, copy_logs



rule prep_multiqc_data:
    input:
        bcbio_mq_filelist = join(run.date_dir, 'multiqc/list_files_final.txt'),
        bcbio_mq_yaml     = join(run.date_dir, 'multiqc/multiqc_config.yaml'),
        bcbio_final_dir   = run.final_dir
    output:
        filelist        = 'work/{batch}/multiqc_data/multiqc_filelist.txt',
        conf_yaml       = 'work/{batch}/multiqc_data/umccr_multiqc_config.yaml',
        bcbio_conf_yaml = 'work/{batch}/multiqc_data/bcbio_multiqc_config.yaml'
    params:
        data_dir        = 'work/{batch}/multiqc_data'
    run:
        conf, additional_files = make_report_metadata(run, base_dirpath=abspath('.'))
        # prep gold standard with:
        # for f in $(grep -v trimmed list_files_final.txt | grep -v target_info.yaml) ; do mkdir -p tmp/`dirname $f` ; cp $f tmp/$f ; done
        gold_standard_dir = join(package_path(), 'multiqc', 'gold_standard')
        with open(join(gold_standard_dir, 'list_files_final.txt')) as f:
            for l in f:
                l = l.strip()
                additional_files.append(join(gold_standard_dir, 'final', l))
        multiqc_prep_data(
            bcbio_mq_filelist=input.bcbio_mq_filelist,
            bcbio_mq_yaml=input.bcbio_mq_yaml,
            bcbio_final_dir=input.bcbio_final_dir,
            new_mq_data_dir=params.data_dir,
            conf=conf,
            filelist_file=output.filelist,
            conf_yaml=output.conf_yaml,
            new_bcbio_mq_yaml=output.bcbio_conf_yaml,
            additional_files=additional_files,
            gold_standard_data=[]
        )


rule batch_multiqc:  # {}
    input:
        filelist        = 'work/{batch}/multiqc_data/multiqc_filelist.txt',
        conf_yaml       = 'work/{batch}/multiqc_data/umccr_multiqc_config.yaml',
        bcbio_conf_yaml = 'work/{batch}/multiqc_data/bcbio_multiqc_config.yaml'
    output:
        html_file       = '{batch}/{batch}-multiqc_report.html'
    run:
        ignore_samples=[s.name for s in run.samples if s.name not in
                [batch_by_name[wildcards.batch].tumor.name, batch_by_name[wildcards.batch].normal.name]]
        ignore_samples_re = '"' + '|'.join(ignore_samples) + '"'
        shell(f'multiqc -f -v -o . -l {input.filelist} -c {input.conf_yaml} -c {input.bcbio_conf_yaml}'
            ' --filename {output.html_file} --ignore-samples {ignore_samples_re}')

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
