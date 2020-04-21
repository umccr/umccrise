from os.path import join

import yaml
from hpc_utils import hpc
from ngs_utils.file_utils import verify_file

localrules: oncoviruses


rule run_oncoviruses:
    input:
        tumor_bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
    output:
        prioritized_tsv = '{batch}/oncoviruses/prioritized_oncoviruses.tsv'
    params:
        genomes_dir = hpc.genomes_dir,
        work_dir = 'work/{batch}/oncoviruses',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
        breakpoints_vcf = '{batch}/oncoviral_breakpoints.vcf.gz',
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
    threads:
        threads_per_batch
    resources:
        mem_mb=10000
    benchmark:
        'benchmarks/{batch}/oncoviruses/{batch}-oncoviruses.tsv'
    group: "oncoviruses"
    run:
        shell('oncoviruses {input.tumor_bam} -o {params.work_dir} -s {params.tumor_name} '
              '--genomes-dir {params.genomes_dir} {params.unlock_opt}')
        if verify_file(join(params.work_dir, 'breakpoints.vcf.gz')):
            shell('cp {params.work_dir}/breakpoints.vcf.gz {params.breakpoints_vcf}')
        shell('cp {params.work_dir}/prioritized_oncoviruses.tsv {output.prioritized_tsv}')


def make_oncoviral_mqc_metric(prioritized_oncoviruses_tsv):
    min_significant_completeness = 0.5
    completeness_threshold = '5x'

    with open(prioritized_oncoviruses_tsv) as f:
        for l in f:
            if l.startswith('##'):
                if l.startswith('## Sample'):
                    sample_name = l.strip().split()[-1]
                continue
            if l.startswith('#'):
                headers = l.strip().split("\t")  # #virus  size    depth   1x      5x      25x
                continue
            viral_data = []
            values_dict = dict(zip(headers, l.strip().split("\t")))
            virus_name = values_dict['#virus']
            ave_depth = float(values_dict['depth'])
            completeness = float(values_dict[completeness_threshold])
            viral_data.append((completeness, ave_depth, virus_name))
    viral_data.sort(reverse=True)
    there_some_hits = any(c >= min_significant_completeness for c, d, v in viral_data)
    if not there_some_hits:
        # showing all that significant, but at least 3 records even if nothing is significant
        viral_data = viral_data[:2]
    else:
        viral_data = [(c, d, v) for c, d, v in viral_data if c >= min_significant_completeness]

    # To add into stats:
    data = {
        sample_name: {
            "viral_content": "; ".join([v for i, (c, d, v) in enumerate(viral_data)]) if there_some_hits else '-'
        }
    }

    header = {
        'title': 'Viral',
        'description': (
            f'Sequences of known oncoviruses, found in umapped reads.'
            f' Format: (x depth; % of sequence covered at >{completeness_threshold}).'
            f' Viral sequences are from'
            f' <a href="https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files">GDC</a>' 
            f' found in unmapped reads. Showing significant hits with at least {completeness_threshold} support' 
            f' along at least {int(100 * min_significant_completeness)}% of the genome.'
        )
    }
    return data, header


rule oncoviral_multiqc:
    input:
        prioritized_tsv = '{batch}/oncoviruses/prioritized_oncoviruses.tsv'
    output:
        yml = 'work/{batch}/oncoviruses/{batch}_oncoviruses_stats.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name
    group: "oncoviruses"
    run:
        data, header = make_oncoviral_mqc_metric(input.prioritized_tsv)

        with open(output.yml, 'w') as out:
            data = {
                'data': data,
                'title': header['title'],
                'description': header['description'],
            }
            yaml.dump(data, out, default_flow_style=False)


rule oncoviruses:
    input:
        expand(rules.run_oncoviruses.output, batch=batch_by_name.keys())
    output:
        temp(touch('log/oncoviruses.done'))
