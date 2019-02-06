from umccrise.multiqc.prep_data import make_report_metadata, multiqc_prep_data
import cyvcf2
import yaml


# localrules: multiqc, copy_logs
localrules: multiqc, copy_config, copy_logs


rule copy_config:
    input:
        conf_dir = run.config_dir
    output:
        conf_dir = directory(join('log/config'))
    shell:
        'cp -r {input.conf_dir} {output.conf_dir}'


# rule copy_logs:
#     input:
#         versions = versions,
#         programs = programs,
#     output:
#         versions = 'log/data_versions.csv',
#         programs = 'log/programs.txt',
#     shell:
#         'cp -r {input.versions} {output.versions} && ' \
#         'cp -r {input.programs} {output.programs}'

rule prep_multiqc_data:
    input:
        bcbio_mq_filelist       = join(run.date_dir, 'multiqc/list_files_final.txt'),
        bcbio_final_dir         = run.final_dir,
        conpair_concord         = rules.run_conpair.output.concord,
        conpair_contam          = rules.run_conpair.output.contam,
        somatic_stats           = rules.somatic_stats_report.output[0],
        germline_stats          = rules.germline_stats_report.output[0],
        bcftools_somatic_stats  = rules.bcftools_stats_somatic.output[0],
    output:
        filelist            = 'work/{batch}/multiqc_data/filelist.txt',
        generated_conf_yaml = 'work/{batch}/multiqc_data/generated_conf.yaml',
    params:
        data_dir        = 'work/{batch}/multiqc_data',
        tumor_name      = lambda wc: batch_by_name[wc.batch].tumor.name,
        normal_name     = lambda wc: batch_by_name[wc.batch].normal.name,
        prog_versions   = join(run.date_dir, 'programs.txt'),
        data_versions   = join(run.date_dir, 'data_versions.csv'),
    # group: 'multiqc'
    run:
        report_base_path = dirname(abspath(f'{wildcards.batch}/{wildcards.batch}-multiqc_report.html'))
        generated_conf, additional_files = make_report_metadata(
            run,
            tumor_sample=params.tumor_name,
            normal_sample=params.normal_name,
            base_dirpath=report_base_path,
            analysis_dir=run.date_dir,
            prog_versions_fpath=verify_file(params.prog_versions, silent=True),
            data_versions_fpath=verify_file(params.data_versions, silent=True),
            new_dir_for_versions=abspath(join(f'{wildcards.batch}', 'log')),
        )
        gold_standard_dir = join(package_path(), 'multiqc', 'gold_standard', 'umccrised_2019.qconly.renamed')
        with open(join(gold_standard_dir, 'background_multiqc_filelist.txt')) as f:
            for l in f:
                l = l.strip()
                additional_files.append(join(gold_standard_dir, l))

        additional_files.extend([
            join(input.conpair_contam, params.tumor_name + '.txt'),
            join(input.conpair_contam, params.normal_name + '.txt'),
            join(input.conpair_concord, params.tumor_name + '.txt'),
            input.somatic_stats,
            input.germline_stats,
            input.bcftools_somatic_stats
        ])

        multiqc_prep_data(
            bcbio_mq_filelist=input.bcbio_mq_filelist,
            bcbio_final_dir=input.bcbio_final_dir,
            new_mq_data_dir=params.data_dir,
            generated_conf=generated_conf,
            out_filelist_file=output.filelist,
            out_conf_yaml=output.generated_conf_yaml,
            additional_files=additional_files,
            exclude_files=[
                '.*qsignature.*',
                '.*bcftools_stats.txt',  # adding somatic stats, but keeping germline stats
                '.*indexcov.tsv',
                '.*ped_check.rel-difference.csv',
                '.*sort-chr.qsig.vcf.*',
                '.*Per_base_N_content.tsv',
                '.*Per_base_sequence_content.tsv',
                '.*Per_base_sequence_quality.tsv',
                '.*Per_sequence_GC_content.tsv',
                '.*Per_sequence_quality_scores.tsv',
                '.*Per_tile_sequence_quality.tsv',
                '.*Sequence_Length_Distribution.tsv',
                '.*.html',
                '.*verifybamid.*',
            ],
            include_files=[
                f'.*{params.tumor_name}.*',
                f'.*{params.normal_name}.*',
                f'.*{wildcards.batch}.*',
            ],
        )


rule batch_multiqc:
    input:
        filelist            = 'work/{batch}/multiqc_data/filelist.txt',
        generated_conf_yaml = 'work/{batch}/multiqc_data/generated_conf.yaml',
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
              f' -c {input.umccrise_conf_yaml} -c {input.generated_conf_yaml} --filename {output.html_file}')


rule multiqc:
    input:
        expand(rules.batch_multiqc.output, batch=batch_by_name.keys()),
        rules.copy_config.output,
    output:
        temp(touch('log/multiqc.done'))













