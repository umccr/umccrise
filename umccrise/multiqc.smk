from umccrise.multiqc.prep_data import make_report_metadata, multiqc_prep_data, parse_bcbio_filelist
import cyvcf2
import yaml
from os.path import abspath, join, dirname
from ngs_utils.file_utils import verify_file
from ngs_utils.bcbio import BcbioProject
from umccrise import package_path


localrules: multiqc


def load_background_samples(genome_build, project_type ='dragen'):
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


if isinstance(run, BcbioProject):
    localrules: copy_config

    # Copy all files from the config/ directory, omitting all directories
    # Omitting directories prevents from errors on IAP when the input mount is read-only,
    # and the directory copied with cp -r inherits the read-only flags, preventing snakemake
    # from creating files inside of it
    rule copy_config:
        input:
            conf_dir = run.config_dir
        output:
            done_flag = 'log/config/.done',
        params:
            conf_dir = 'log/config',
        shell:
            'for f in {input.conf_dir}/*; do test ! -f $f || cp $f {params.conf_dir}/; done && '
            'touch {output.done_flag}'

    rule prep_multiqc_data:
        input:
            bcbio_mq_filelist       = join(run.date_dir, 'multiqc/list_files_final.txt'),
            bcbio_final_dir         = run.final_dir,
            conpair_concord         = rules.run_conpair.output.concord,
            conpair_contam          = rules.run_conpair.output.contam,
            somatic_stats           = rules.somatic_stats_report.output[0],
            germline_stats          = rules.germline_stats_report.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
            bcftools_somatic_stats  = rules.bcftools_stats_somatic.output[0],
            bcftools_germline_stats = rules.bcftools_stats_germline.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
            prog_versions           = join(run.date_dir, 'programs.txt'),
            data_versions           = join(run.date_dir, 'data_versions.csv'),
            # oncoviruses_list        = rules.run_oncoviruses.output.prioritized_tsv,
            oncoviruses_stats       = rules.oncoviral_multiqc.output.yml,
        output:
            filelist                = 'work/{batch}/multiqc_data/filelist.txt',
            generated_conf_yaml     = 'work/{batch}/multiqc_data/generated_conf.yaml',
            prog_versions           = '{batch}/log/programs.txt',
            data_versions           = '{batch}/log/data_versions.txt',
        params:
            data_dir                = 'work/{batch}/multiqc_data',
            genome_build            = run.genome_build
        group: 'multiqc'
        run:
            batch = batch_by_name[wildcards.batch]
            report_base_path = dirname(abspath(f'{wildcards.batch}/{wildcards.batch}-multiqc_report.html'))

            generated_conf, qc_files = make_report_metadata(
                run,
                batch=batch,
                base_dirpath=report_base_path,
                analysis_dir=run.date_dir,
                prog_versions_fpath=verify_file(input.prog_versions, silent=True),
                data_versions_fpath=verify_file(input.data_versions, silent=True),
                new_dir_for_versions=dirname(output.prog_versions)
            )

            # Gold standard QC files
            qc_files.extend(load_background_samples(params.genome_build, project_type ='bcbio'))

            # Umccrise QC files
            qc_files.extend([
                join(input.conpair_concord, batch.tumor.name + '.concordance.txt'),
                join(input.conpair_contam, batch.normal.name + '.contamination.txt'),
                join(input.conpair_contam, batch.tumor.name + '.contamination.txt'),
                input.somatic_stats,
                input.germline_stats,
                input.bcftools_somatic_stats,
                input.bcftools_germline_stats,
                input.oncoviruses_stats,
            ])

            # Bcbio QC files
            qc_files.extend(parse_bcbio_filelist(
                bcbio_mq_filelist=input.bcbio_mq_filelist,
                bcbio_final_dir=input.bcbio_final_dir,
                new_mq_data_dir=params.data_dir,
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
                    '.*viral.*',
                ],
                include_files=[
                    f'.*{batch.tumor.name}.*',
                    f'.*{batch.normal.name}.*',
                    f'.*{batch.name}.*',
                ],
            ))

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
            other_samples=[
                s.name for s in run.samples if s.name not in [batch_by_name[wildcards.batch].tumor.name,
                                                              batch_by_name[wildcards.batch].normal.name]]
            if other_samples:
                greps = ''.join(f' | grep -v "__{sn}/" | grep -v "/{sn}/" | grep -v "/{sn}_bcbio.txt"' for sn in other_samples)
                list_files = f'<(cat {input.filelist}{greps})'
            else:
                list_files = input.filelist
            shell(f'LC_ALL=$LC_ALL LANG=$LANG multiqc -f -o . -l {list_files}'
                  f' -c {input.umccrise_conf_yaml} -c {input.generated_conf_yaml} --filename {output.html_file}')


else:  # dragen
    rule prep_multiqc_data:
        input:
            conpair_concord         = rules.run_conpair.output.concord,
            conpair_contam          = rules.run_conpair.output.contam,
            somatic_stats           = rules.somatic_stats_report.output[0],
            germline_stats          = rules.germline_stats_report.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
            bcftools_somatic_stats  = rules.bcftools_stats_somatic.output[0],
            bcftools_germline_stats = rules.bcftools_stats_germline.output[0] if all(b.germline_vcf for b in batch_by_name.values()) else [],
        output:
            filelist                = 'work/{batch}/multiqc_data/filelist.txt',
            generated_conf_yaml     = 'work/{batch}/multiqc_data/generated_conf.yaml',
        params:
            data_dir                = 'work/{batch}/multiqc_data',
            genome_build            = run.genome_build
        group: 'multiqc'
        run:
            batch = batch_by_name[wildcards.batch]
            report_base_path = dirname(abspath(f'{wildcards.batch}/{wildcards.batch}-multiqc_report.html'))
            generated_conf, qc_files = make_report_metadata(
                run,
                batch=batch,
                base_dirpath=report_base_path,
            )

            # Gold standard QC files
            qc_files.extend(load_background_samples(params.genome_build, project_type ='dragen'))

            # Umccrise QC files
            qc_files.extend([
                join(input.conpair_concord, batch.tumor.name + '.concordance.txt'),
                join(input.conpair_contam, batch.normal.name + '.contamination.txt'),
                join(input.conpair_contam, batch.tumor.name + '.contamination.txt'),
                input.somatic_stats,
                input.germline_stats,
                input.bcftools_somatic_stats,
                input.bcftools_germline_stats,
            ])

            # Dragen QC files
            qc_files.extend(batch_by_name[wildcards.batch].all_qc_files()),

            multiqc_prep_data(
                generated_conf=generated_conf,
                out_filelist_file=output.filelist,
                out_conf_yaml=output.generated_conf_yaml,
                qc_files=qc_files,
            )

    rule batch_multiqc:
        input:
            umccrise_conf_yaml  = join(package_path(), 'multiqc', 'multiqc_config.yaml'),
            filelist            = 'work/{batch}/multiqc_data/filelist.txt',
            generated_conf_yaml = 'work/{batch}/multiqc_data/generated_conf.yaml',
        output:
            html_file = '{batch}/{batch}-multiqc_report.html'
        group: 'multiqc'
        run:
            list_files = input.filelist
            shell(f'LC_ALL=$LC_ALL LANG=$LANG multiqc -f -o . -l {list_files}'
                  f' -c {input.umccrise_conf_yaml} -c {input.generated_conf_yaml} --filename {output.html_file}')


if len(batch_by_name) > 1:
    rule combined_multiqc_prep_multiqc_data:
        input:
            filelists           = [f'work/{b}/multiqc_data/filelist.txt' for b in batch_by_name.keys()],
        output:
            filelist            = 'multiqc/multiqc_data/filelist.txt',
            generated_conf_yaml = 'multiqc/multiqc_data/generated_conf.yaml',
        group: 'combined_multiqc'
        run:
            report_base_path = dirname(abspath(f'multiqc/multiqc_report.html'))
            generated_conf, _ = make_report_metadata(
                run,
                base_dirpath=report_base_path,
            )

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
        rules.copy_config.output if isinstance(run, BcbioProject) else [],
        rules.combined_multiqc.output if len(batch_by_name) > 1 else [],
    output:
        temp(touch('log/multiqc.done'))









