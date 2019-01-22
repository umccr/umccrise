from umccrise.prep_multiqc import make_report_metadata, multiqc_prep_data


# localrules: multiqc, copy_logs
localrules: prep_multiqc_data, multiqc, copy_logs


"""
#############################
### Prepare gold standard ###
#############################
cd /data/cephfs/punim0010/data/Results/Tothill-A5/2018-08-11
mkdir final.subset
for f in $(cat final/2018-08-11_2018-07-31T0005_Tothill-A5_WGS-merged/multiqc/list_files_final.txt | grep -v trimmed | grep -v target_info.yaml) ; do 
    mkdir -p final.subset/`dirname $f`
    cp final/$f final.subset/$f 
done

cat final/2018-08-11_2018-07-31T0005_Tothill-A5_WGS-merged/multiqc/list_files_final.txt |\
    grep -P "E201|E199|E194|E190|E202" |\
    grep -v "indexcov.tsv" |\
    grep -v "Per_base_N_content.tsv" |\
    grep -v "Per_base_sequence_content.tsv" |\
    grep -v "Per_base_sequence_quality.tsv" |\
    grep -v "Per_sequence_GC_content.tsv" |\
    grep -v "Per_sequence_quality_scores.tsv" |\
    grep -v "Per_tile_sequence_quality.tsv" |\
    grep -v "Sequence_Length_Distribution.tsv" |\
    grep -v "sort-chr.qsig.vcf" |\
    grep -v "ped_check.rel-difference.csv" |\
    grep -v ".html" |\
    > final.subset/list_files_final.txt

### Clean up ###
cd final.subset

# Remove sample dirs:
ls | grep -v -P "E201|E199|E194|E190|E202|2018-08-11_2018-07-31T0005_Tothill-A5_WGS-merged|list_files_final.txt" | xargs rm -rf

# Clean up qc/coverage
find . -name "*-indexcov.tsv" -delete
# Clean up qc/fastqc
find . -name "fastqc_report.html" -delete
find . -name "Per_base_N_content.tsv" -delete
find . -name "Per_base_sequence_content.tsv" -delete
find . -name "Per_base_sequence_quality.tsv" -delete
find . -name "Per_sequence_GC_content.tsv" -delete
find . -name "Per_sequence_quality_scores.tsv" -delete
find . -name "Per_tile_sequence_quality.tsv" -delete
find . -name "Sequence_Length_Distribution.tsv" -delete
# Clean up qc/qsignature
rm -rf */qc/qsignature
# Clean up qc/peddy -  TODO: remove peddy and VerifyBAMID
find . -path "*qc/peddy/*.ped_check.rel-difference.csv" -delete
find . -path "*qc/peddy/*.html" -delete
# Clean up bcbio metrics
cd 2018-08-11_2018-07-31T0005_Tothill-A5_WGS-merged/multiqc/report/metrics
ls | grep -v -P "E201|E199|E194|E190|E202" | xargs rm 
cd ../../../../

cd ..
### Rename ###
python rename.py final.subset final.subset.renamed
"""


def prep_conpair(contam, concor, data_dir, t_name, n_name):
    contam_dir = safe_mkdir(join(data_dir, 'conpair', 'contamination'))
    concor_dir = safe_mkdir(join(data_dir, 'conpair', 'concordance'))
    contam_t = join(contam_dir, t_name + '.txt')
    contam_n = join(contam_dir, n_name + '.txt')

    with open(contam) as c, open(contam_t, 'w') as ct, open(contam_n, 'w') as cn:
        for l in c:
            if l.startswith('Tumor sample contamination level'):
                ct.write(l.replace('Tumor sample', 'Sample'))
            if l.startswith('Normal sample contamination level'):
                cn.write(l.replace('Normal sample', 'Sample'))

    concor_t = join(concor_dir, t_name + '.txt')
    shell(f"cp {concor} {concor_t}")

    return contam_t, contam_n, concor_t


rule prep_multiqc_data:
    input:
        bcbio_mq_filelist = join(run.date_dir, 'multiqc/list_files_final.txt'),
        bcbio_mq_yaml     = join(run.date_dir, 'multiqc/multiqc_config.yaml'),
        bcbio_final_dir   = run.final_dir,
        versions          = 'log/' + run.project_name + '-data_versions.csv',
        programs          = 'log/' + run.project_name + '-programs.txt',
        conpair_contam    = '{batch}/conpair/tumor_normal_contamination.txt',
        conpair_concor    = '{batch}/conpair/tumor_normal_concordance.txt',
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

        additional_files.extend(prep_conpair(input.conpair_contam, input.conpair_concor,
                                             params.data_dir, params.tumor_name, params.normal_name))

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












