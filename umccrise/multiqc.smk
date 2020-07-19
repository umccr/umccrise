import glob
from ngs_utils.utils import update_dict
from umccrise.multiqc.prep_data import make_report_metadata, multiqc_prep_data
import yaml
from os.path import abspath, join, dirname, basename
from ngs_utils.file_utils import verify_file
from ngs_utils.bcbio import BcbioProject, BcbioBatch
from ngs_utils.dragen import DragenProject
from umccrise import package_path


localrules: multiqc


def load_background_samples(genome_build, project_type='dragen'):
    paths = []
    if project_type == 'bcbio':
        gold_standard_dir = join(package_path(), 'multiqc', 'gold_standard', 'umccrised.qconly.renamed')
    else:
        gold_standard_dir = join(package_path(), 'multiqc', 'gold_standard', 'dragen', 'umccrised.qconly')
    if genome_build == 'hg38':
        gold_standard_dir += '_hg38'
    with open(join(gold_standard_dir, 'background_multiqc_filelist.txt')) as f:
        for l in f:
            l = l.strip()
            if l:
                paths.append(join(gold_standard_dir, l))
    return paths


def find_bcbio_qc_files(batch: BcbioBatch, dst_dir):
    """
    Finds all QC files by parsing the MultiQC file_list, copies them into `dst_dir`
    :param bcbio_run: BcbioProject object
    :param dst_dir: destination directory where the QC files will be copied to
    :return: list of found QC files
    """
    bcbio_qc_files = batch.find_qc_files(
        dst_dir,
        # exclude and include applies only to the sample files, not to the background.
        # first excluding, then including from the remaining.
        exclude_files=[
            '.*qsignature.*',
            '.*bcftools_stats*.txt',
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
            '.*gdc-viral-completeness.txt',
        ] + (['.*peddy.*'] if 'peddy' in stages else []),
        include_files=[
            f'.*{batch.name}.*',
            f'.*{batch.tumor.name}.*',
        ] + [f'.*{batch.normal.name}.*'] if batch.normal else [],
    )
    return bcbio_qc_files


rule prep_multiqc_data:
    input:
        conpair_concord         = rules.run_conpair.output.concord,
        conpair_contam          = rules.run_conpair.output.contam,
        somatic_stats           = rules.somatic_stats_report.output[0],
        germline_stats          = rules.germline_stats_report.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
        bcftools_somatic_stats  = rules.bcftools_stats_somatic.output[0],
        bcftools_germline_stats = rules.bcftools_stats_germline.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
        oncoviruses_data        = rules.oncoviral_multiqc.output.data_yml,
        oncoviruses_header      = rules.oncoviral_multiqc.output.header_yml,
        purple_stats            = rules.purple_run.output.purity,
        purple_qc               = rules.purple_run.output.qc,
        peddy_dirs              = expand(rules.run_peddy.output.dir.replace('{batch}', '{{batch}}'), phenotype=['normal'])
                                    if 'peddy' in stages else [],
    output:
        filelist                = 'work/{batch}/multiqc_data/filelist.txt',
        generated_conf_yaml     = 'work/{batch}/multiqc_data/generated_conf.yaml',
    params:
        data_dir                = 'work/{batch}/multiqc_data',
        genome_build            = run.genome_build,
    group: 'multiqc'
    run:
        batch = batch_by_name[wildcards.batch]
        report_base_path = dirname(abspath(f'{wildcards.batch}/{wildcards.batch}-multiqc_report.html'))

        if isinstance(run, BcbioProject):
            prog_versions     = join(run.date_dir, 'programs.txt')
            data_versions     = join(run.date_dir, 'data_versions.csv')

            generated_conf, qc_files = make_report_metadata(
                run,
                batch=batch,
                base_dirpath=report_base_path,
                prog_versions_fpath=verify_file(prog_versions, silent=True),
                data_versions_fpath=verify_file(data_versions, silent=True),
                new_dir_for_versions=f'{wildcards.batch}/log'
            )
        else:
            generated_conf, qc_files = make_report_metadata(
                run,
                batch=batch,
                base_dirpath=report_base_path,
            )

        with open(input.oncoviruses_header) as f:
            generated_conf.update(yaml.load(f))

        # Gold standard QC files
        if isinstance(batch, BcbioBatch):
            qc_files.extend(load_background_samples(params.genome_build, project_type='bcbio'))
        if isinstance(batch, DragenProject):
            qc_files.extend(load_background_samples(params.genome_build, project_type='dragen'))

        # Umccrise QC files
        renamed_purple_qc    = join(params.data_dir, basename(input.purple_qc   ).replace(wildcards.batch, batch.tumor.name))
        renamed_purple_stats = join(params.data_dir, basename(input.purple_stats).replace(wildcards.batch, batch.tumor.name))
        shell(f'cp {input.purple_qc} {renamed_purple_qc}')
        shell(f'cp {input.purple_stats} {renamed_purple_stats}')
        qc_files.extend([
            join(input.conpair_concord, batch.tumor.name + '.concordance.txt'),
            join(input.conpair_contam, batch.normal.name + '.contamination.txt'),
            join(input.conpair_contam, batch.tumor.name + '.contamination.txt'),
            input.somatic_stats,
            input.germline_stats,
            input.bcftools_somatic_stats,
            input.bcftools_germline_stats,
            input.oncoviruses_data,
            renamed_purple_qc,
            renamed_purple_stats,
        ])
        peddy_files = []
        for d in input.peddy_dirs:
            peddy_files.extend(glob.glob(join(d, '*.peddy.ped')))
            peddy_files.extend(glob.glob(join(d, '*.het_check.csv')))
            peddy_files.extend(glob.glob(join(d, '*.ped_check.csv')))
            peddy_files.extend(glob.glob(join(d, '*.sex_check.csv')))
        qc_files.extend(peddy_files)

        if isinstance(batch, BcbioBatch):
            qc_files.extend(find_bcbio_qc_files(batch, params.data_dir))
        elif isinstance(run, DragenProject):
            qc_files.extend(batch.all_qc_files()),

        multiqc_prep_data(
            generated_conf=generated_conf,
            out_filelist_file=output.filelist,
            out_conf_yaml=output.generated_conf_yaml,
            qc_files=qc_files
        )


rule batch_multiqc:
    input:
        filelist            = 'work/{batch}/multiqc_data/filelist.txt',
        generated_conf_yaml = 'work/{batch}/multiqc_data/generated_conf.yaml',
        umccrise_conf_yaml  = join(package_path(), 'multiqc', 'multiqc_config.yaml'),
    output:
        html_file           = '{batch}/{batch}-multiqc_report.html'
    group: 'multiqc'
    run:
        list_files = input.filelist
        if isinstance(run, BcbioProject):
            other_samples=[
                s.name for s in run.samples if s.name not in [batch_by_name[wildcards.batch].tumor.name,
                                                              batch_by_name[wildcards.batch].normal.name]]
            if other_samples:
                greps = ''.join(f' | grep -v "__{sn}/" | grep -v "/{sn}/" | grep -v "/{sn}_bcbio.txt"' for sn in other_samples)
                list_files = f'<(cat {input.filelist}{greps})'
        shell(f'LC_ALL=$LC_ALL LANG=$LANG multiqc -f -o . -l {list_files}'
              f' -c {input.umccrise_conf_yaml} -c {input.generated_conf_yaml} --filename {output.html_file}')


# if isinstance(run, BcbioProject):
#     rule prep_multiqc_data:
#         input:
#             bcbio_mq_filelist       = join(run.date_dir, 'multiqc/list_files_final.txt'),
#             bcbio_final_dir         = run.final_dir,
#             conpair_concord         = rules.run_conpair.output.concord,
#             conpair_contam          = rules.run_conpair.output.contam,
#             somatic_stats           = rules.somatic_stats_report.output[0],
#             germline_stats          = rules.germline_stats_report.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
#             bcftools_somatic_stats  = rules.bcftools_stats_somatic.output[0],
#             bcftools_germline_stats = rules.bcftools_stats_germline.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
#             prog_versions           = join(run.date_dir, 'programs.txt'),
#             data_versions           = join(run.date_dir, 'data_versions.csv'),
#             oncoviruses_data        = rules.oncoviral_multiqc.output.data_yml,
#             oncoviruses_header      = rules.oncoviral_multiqc.output.header_yml,
#             purple_stats            = rules.purple_run.output.purity,
#             purple_qc               = rules.purple_run.output.qc,
#         output:
#             filelist                = 'work/{batch}/multiqc_data/filelist.txt',
#             generated_conf_yaml     = 'work/{batch}/multiqc_data/generated_conf.yaml',
#             prog_versions           = '{batch}/log/programs.txt',
#             data_versions           = '{batch}/log/data_versions.txt',
#         params:
#             data_dir                = 'work/{batch}/multiqc_data',
#             genome_build            = run.genome_build
#         group: 'multiqc'
#         run:
#             batch = batch_by_name[wildcards.batch]
#             report_base_path = dirname(abspath(f'{wildcards.batch}/{wildcards.batch}-multiqc_report.html'))
#
#             generated_conf, qc_files = make_report_metadata(
#                 run,
#                 batch=batch,
#                 base_dirpath=report_base_path,
#                 analysis_dir=run.date_dir,
#                 prog_versions_fpath=verify_file(input.prog_versions, silent=True),
#                 data_versions_fpath=verify_file(input.data_versions, silent=True),
#                 new_dir_for_versions=dirname(output.prog_versions)
#             )
#             with open(input.oncoviruses_header) as f:
#                 generated_conf.update(yaml.load(f))
#
#             # Gold standard QC files
#             qc_files.extend(load_background_samples(params.genome_build, project_type ='bcbio'))
#
#             renamed_purple_qc    = join(params.data_dir, basename(input.purple_qc   ).replace(wildcards.batch, batch.tumor.name))
#             renamed_purple_stats = join(params.data_dir, basename(input.purple_stats).replace(wildcards.batch, batch.tumor.name))
#             shell(f'cp {input.purple_qc} {renamed_purple_qc}')
#             shell(f'cp {input.purple_stats} {renamed_purple_stats}')
#
#             # Umccrise QC files
#             qc_files.extend([
#                 join(input.conpair_concord, batch.tumor.name + '.concordance.txt'),
#                 join(input.conpair_contam, batch.normal.name + '.contamination.txt'),
#                 join(input.conpair_contam, batch.tumor.name + '.contamination.txt'),
#                 input.somatic_stats,
#                 input.germline_stats,
#                 input.bcftools_somatic_stats,
#                 input.bcftools_germline_stats,
#                 input.oncoviruses_data,
#                 renamed_purple_qc,
#                 renamed_purple_stats,
#             ])
#
#             # Bcbio QC files
#             bcbio_qc_files = parse_bcbio_filelist(
#                 bcbio_mq_filelist=input.bcbio_mq_filelist,
#                 bcbio_final_dir=input.bcbio_final_dir,
#                 new_mq_data_dir=params.data_dir,
#                 # exclude and include applies only to the sample files, not to the background.
#                 # first excluding, then including from the remaining.
#                 exclude_files=[
#                     '.*qsignature.*',
#                     '.*bcftools_stats*.txt',
#                     '.*indexcov.tsv',
#                     '.*ped_check.rel-difference.csv',
#                     '.*sort-chr.qsig.vcf.*',
#                     '.*Per_base_N_content.tsv',
#                     '.*Per_base_sequence_content.tsv',
#                     '.*Per_base_sequence_quality.tsv',
#                     '.*Per_sequence_GC_content.tsv',
#                     '.*Per_sequence_quality_scores.tsv',
#                     '.*Per_tile_sequence_quality.tsv',
#                     '.*Sequence_Length_Distribution.tsv',
#                     '.*.html',
#                     '.*verifybamid.*',
#                     '.*gdc-viral-completeness.txt',
#                 ],
#                 include_files=[
#                     f'.*{batch.tumor.name}.*',
#                     f'.*{batch.normal.name}.*',
#                     f'.*{batch.name}.*',
#                 ],
#             )
#             qc_files.extend(bcbio_qc_files)
#
#             multiqc_prep_data(
#                 generated_conf=generated_conf,
#                 out_filelist_file=output.filelist,
#                 out_conf_yaml=output.generated_conf_yaml,
#                 qc_files=qc_files
#             )
#
#
#     rule batch_multiqc:
#         input:
#             filelist            = 'work/{batch}/multiqc_data/filelist.txt',
#             generated_conf_yaml = 'work/{batch}/multiqc_data/generated_conf.yaml',
#             umccrise_conf_yaml  = join(package_path(), 'multiqc', 'multiqc_config.yaml'),
#         output:
#             html_file           = '{batch}/{batch}-multiqc_report.html'
#         group: 'multiqc'
#         run:
#             other_samples=[
#                 s.name for s in run.samples if s.name not in [batch_by_name[wildcards.batch].tumor.name,
#                                                               batch_by_name[wildcards.batch].normal.name]]
#             if other_samples:
#                 greps = ''.join(f' | grep -v "__{sn}/" | grep -v "/{sn}/" | grep -v "/{sn}_bcbio.txt"' for sn in other_samples)
#                 list_files = f'<(cat {input.filelist}{greps})'
#             else:
#                 list_files = input.filelist
#             shell(f'LC_ALL=$LC_ALL LANG=$LANG multiqc -f -o . -l {list_files}'
#                   f' -c {input.umccrise_conf_yaml} -c {input.generated_conf_yaml} --filename {output.html_file}')
#
#
# else:  # dragen or any custom data
#     rule prep_multiqc_data:
#         input:
#             conpair_concord         = rules.run_conpair.output.concord,
#             conpair_contam          = rules.run_conpair.output.contam,
#             somatic_stats           = rules.somatic_stats_report.output[0],
#             germline_stats          = rules.germline_stats_report.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
#             bcftools_somatic_stats  = rules.bcftools_stats_somatic.output[0],
#             bcftools_germline_stats = rules.bcftools_stats_germline.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
#             oncoviruses_data        = rules.oncoviral_multiqc.output.data_yml,
#             oncoviruses_header      = rules.oncoviral_multiqc.output.header_yml,
#             purple_stats            = rules.purple_run.output.purity,
#             purple_qc               = rules.purple_run.output.qc,
#         output:
#             filelist                = 'work/{batch}/multiqc_data/filelist.txt',
#             generated_conf_yaml     = 'work/{batch}/multiqc_data/generated_conf.yaml',
#         params:
#             data_dir                = 'work/{batch}/multiqc_data',
#             genome_build            = run.genome_build
#         group: 'multiqc'
#         run:
#             batch = batch_by_name[wildcards.batch]
#             report_base_path = dirname(abspath(f'{wildcards.batch}/{wildcards.batch}-multiqc_report.html'))
#
#             generated_conf, qc_files = make_report_metadata(
#                 run,
#                 batch=batch,
#                 base_dirpath=report_base_path,
#             )
#             with open(input.oncoviruses_header) as f:
#                 generated_conf.update(yaml.load(f))
#
#             # Gold standard QC files
#             qc_files.extend(load_background_samples(params.genome_build, project_type ='dragen'))
#
#             renamed_purple_qc    = join(params.data_dir, basename(input.purple_qc   ).replace(batch.name, batch.tumor.name))
#             renamed_purple_stats = join(params.data_dir, basename(input.purple_stats).replace(batch.name, batch.tumor.name))
#             shell(f'cp {input.purple_qc} {renamed_purple_qc}')
#             shell(f'cp {input.purple_stats} {renamed_purple_stats}')
#
#             # Umccrise QC files
#             qc_files.extend([
#                 join(input.conpair_concord, batch.tumor.name + '.concordance.txt'),
#                 join(input.conpair_contam, batch.normal.name + '.contamination.txt'),
#                 join(input.conpair_contam, batch.tumor.name + '.contamination.txt'),
#                 input.somatic_stats,
#                 input.germline_stats,
#                 input.bcftools_somatic_stats,
#                 input.bcftools_germline_stats,
#                 input.oncoviruses_data,
#                 renamed_purple_qc,
#                 renamed_purple_stats,
#             ])
#
#             # Dragen QC files
#             qc_files.extend(batch_by_name[wildcards.batch].all_qc_files()),
#
#             multiqc_prep_data(
#                 generated_conf=generated_conf,
#                 out_filelist_file=output.filelist,
#                 out_conf_yaml=output.generated_conf_yaml,
#                 qc_files=qc_files,
#             )
#
#     rule batch_multiqc:
#         input:
#             umccrise_conf_yaml  = join(package_path(), 'multiqc', 'multiqc_config.yaml'),
#             filelist            = 'work/{batch}/multiqc_data/filelist.txt',
#             generated_conf_yaml = 'work/{batch}/multiqc_data/generated_conf.yaml',
#         output:
#             html_file = '{batch}/{batch}-multiqc_report.html'
#         group: 'multiqc'
#         run:
#             list_files = input.filelist
#             shell(f'LC_ALL=$LC_ALL LANG=$LANG multiqc -f -o . -l {list_files}'
#                   f' -c {input.umccrise_conf_yaml} -c {input.generated_conf_yaml} --filename {output.html_file}')


if len(batch_by_name) > 1:
    rule combined_multiqc_prep_multiqc_data:
        input:
            filelists            = [f'work/{b}/multiqc_data/filelist.txt' for b in batch_by_name.keys()],
            generated_conf_yamls = [f'work/{b}/multiqc_data/generated_conf.yaml' for b in batch_by_name.keys()],
        output:
            filelist             = 'multiqc/multiqc_data/filelist.txt',
            generated_conf_yaml  = 'multiqc/multiqc_data/generated_conf.yaml',
        group: 'combined_multiqc'
        run:
            report_base_path = dirname(abspath(f'multiqc/multiqc_report.html'))
            generated_conf, _ = make_report_metadata(
                run,
                base_dirpath=report_base_path
            )
            for sample_conf_yaml in input.generated_conf_yamls:
                with open(sample_conf_yaml) as f:
                    sample_conf = yaml.load(f)
                    del sample_conf['umccr']
                    generated_conf = update_dict(generated_conf, sample_conf)

            qc_files = []
            for filelist in input.filelists:
                with open(filelist) as f:
                    for l in f:
                        l = l.strip()
                        if l and l not in qc_files:
                            qc_files.append(l)

            multiqc_prep_data(
                generated_conf=generated_conf,
                out_filelist_file=output.filelist,
                out_conf_yaml=output.generated_conf_yaml,
                qc_files=qc_files,
            )

    rule combined_multiqc:
        input:
            umccrise_conf_yaml  = join(package_path(), 'multiqc', 'multiqc_config.yaml'),
            filelist            = 'multiqc/multiqc_data/filelist.txt',
            generated_conf_yaml = 'multiqc/multiqc_data/generated_conf.yaml',
        output:
            html_file           = 'multiqc/multiqc_report.html'
        group: 'combined_multiqc'
        run:
            list_files = input.filelist
            shell(f'LC_ALL=$LC_ALL LANG=$LANG multiqc -f -o . -l {list_files}'
                  f' -c {input.umccrise_conf_yaml} -c {input.generated_conf_yaml} --filename {output.html_file}')


rule multiqc:
    input:
        expand(rules.batch_multiqc.output, batch=batch_by_name.keys()),
        rules.combined_multiqc.output if len(batch_by_name) > 1 else [],
    output:
        temp(touch('log/multiqc.done'))


