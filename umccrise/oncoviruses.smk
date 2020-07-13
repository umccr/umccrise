import re
from collections import defaultdict
from os.path import join, basename
import yaml
from cyvcf2 import VCF
import csv
from ngs_utils.file_utils import verify_file
from ngs_utils.utils import update_dict
from reference_data import api as refdata


localrules: oncoviruses, oncoviruses_per_batch


checkpoint viral_content:
    input:
        tumor_bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
    output:
        prioritized_tsv = '{batch}/oncoviruses/prioritized_oncoviruses.tsv',
        present_viruses = 'work/{batch}/oncoviruses/present_viruses.txt',
    params:
        genomes_dir = refdata.genomes_dir,
        work_dir = 'work/{batch}/oncoviruses',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
    threads:
        threads_per_batch
    resources:
        mem_mb = lambda wildcards, attempt: 20000 * attempt,
    benchmark:
        'benchmarks/{batch}/oncoviruses/{batch}-oncoviruses.tsv'
    group: "viral_content"
    run:
        shell('oncoviruses {input.tumor_bam} -o {params.work_dir} -s {params.tumor_name} '
              '--genomes-dir {params.genomes_dir} {params.unlock_opt} --only-detect')
        # if verify_file(join(params.work_dir, 'breakpoints.vcf.gz')):
        #     shell('cp {params.work_dir}/breakpoints.vcf.gz {params.breakpoints_vcf}')
        shell('cp {params.work_dir}/prioritized_oncoviruses.tsv {output.prioritized_tsv}')

        viruses = []
        with open(output.prioritized_tsv) as f:
            for l in f:
                l = l.strip()
                if l and not l.startswith('#'):
                    l = l.split('\t')
                    if l[6] != '.':
                        viruses.append(l[0])
        with open(output.present_viruses, 'w') as outf:
            outf.write(','.join(viruses))


rule viral_integration_sites:
    input:
        tumor_bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
        significant_viruses = 'work/{batch}/oncoviruses/present_viruses.txt',
    output:
        breakpoints_vcf = '{batch}/oncoviruses/oncoviral_breakpoints.vcf.gz',
    params:
        genomes_dir = refdata.genomes_dir,
        work_dir = 'work/{batch}/oncoviruses',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
        tumor_name = lambda wc: batch_by_name[wc.batch].tumor.name,
    threads:
        threads_per_batch
    resources:
        mem_mb = lambda wildcards, attempt: 10000 * attempt,
    benchmark:
        'benchmarks/{batch}/oncoviruses/{batch}-oncoviruses.tsv'
    group: "viral_is"
    shell:
        'oncoviruses {input.tumor_bam} -o {params.work_dir} -s {params.tumor_name} '
        '--genomes-dir {params.genomes_dir} {params.unlock_opt} -v $(cat {input.significant_viruses})'
        '; cp {params.work_dir}/breakpoints.vcf.gz {output.breakpoints_vcf}'


def parse_info_field(rec, name):
    val = rec.INFO.get(name)
    if val is None:
        return ''
    elif isinstance(val, float) or isinstance(val, int) or isinstance(val, bool) or isinstance(val, str):
        return str(val)
    else:
        return ','.join(map(str, val))

# Produce a TSV file for further analysis in Rmd and MultiQC
# sample                chrom  start      end  svtype  PAIR_COUNT  ViralGenes  GenesWithin100kb                               ID                      MATEID
# PRJ180253_E190-T01-D  chr8   127719201  .    BND     111         L1          AC104370.1,AC108925.1,CASC11,MYC,MIR1204,PVT1  MantaBND:0:1:5:0:0:0:1  MantaBND:0:1:5:0:0:0:0
# PRJ180253_E190-T01-D  HPV18       6787  .    BND     111         L1          AC104370.1,AC108925.1,CASC11,MYC,MIR1204,PVT1  MantaBND:0:1:5:0:0:0:0  MantaBND:0:1:5:0:0:0:1
# OR:
# sample                virus  chrom  host_start  virus_start  PAIR_COUNT  ViralGenes  GenesWithin100kb                               ID                      MATEID
# PRJ180253_E190-T01-D  HPV18  chr8   127719201   6787         111         L1          AC104370.1,AC108925.1,CASC11,MYC,MIR1204,PVT1  MantaBND:0:1:5:0:0:0:1  MantaBND:0:1:5:0:0:0:0
rule oncoviruses_breakpoints_tsv:
    input:
        vcf = '{batch}/oncoviruses/oncoviral_breakpoints.vcf.gz',
        present_viruses = 'work/{batch}/oncoviruses/present_viruses.txt',
    output:
        tsv = 'work/{batch}/oncoviruses/oncoviral_breakpoints.tsv'
    group: "viral_is"
    run:
        sample_name = batch_by_name[wildcards.batch].tumor.name
        viruses = open(input.present_viruses).read().split(',')
        with open(output.tsv, 'w') as out:
            header = ['sample', 'contig', 'start', 'end', 'svtype',
                      'PAIR_COUNT', 'DisruptedGenes', 'UpstreamGenes', 'ID', 'MATEID',]
            out.write('\t'.join(header) + '\n')
            for rec in VCF(input.vcf):
                viral_genes = rec.INFO.get('ViralGenes', '').split(',')
                disrupted_host_genes = rec.INFO.get('DisruptedGenes', '').split(',')
                _genes_within_100kb = rec.INFO.get('GenesWithin100kb', '').split(',')
                upstream_host_genes = [g for g in _genes_within_100kb if g not in disrupted_host_genes]
                if 'PE' in rec.FORMAT or 'SR' in rec.FORMAT:
                    read_support = 0
                    if 'PE' in rec.FORMAT:
                        read_support += int(rec.format('PE')[0])
                    if 'SR' in rec.FORMAT:
                        read_support += int(rec.format('SR')[0])
                elif 'PAIR_COUNT' in rec.INFO:
                    read_support = parse_info_field(rec, 'PAIR_COUNT')
                else:
                    read_support = ''
                disrupted_genes = viral_genes if rec.CHROM in viruses else disrupted_host_genes
                data = [
                    sample_name, rec.CHROM, rec.POS, rec.INFO.get('END', ''),
                    rec.INFO['SVTYPE'],
                    str(read_support),
                    ', '.join(disrupted_genes) if disrupted_genes else 'None',
                    '.' if rec.CHROM in viruses else ', '.join(upstream_host_genes),
                    rec.ID,
                    parse_info_field(rec, 'MATEID'),
                ]
                out.write('\t'.join(map(str, data)) + '\n')


def make_oncoviral_mqc_metric(prioritized_oncoviruses_tsv):
    viral_data = []
    with open(prioritized_oncoviruses_tsv) as f:
        ## Viral sequences (from {ONCOVIRAL_SOURCE_URL}) found in unmapped reads
        ## Sample: {SAMPLE}
        ## Minimal completeness: {MIN_1x_PCT}% at 1x or {MIN_5x_LEN}bp at 5x
        #virus\tsize\tdepth\t1x\t5x\t25x\tsignificance
        sname, min_1x_pct, min_5x_len = None, None, None
        for l in f:
            if l.startswith('##'):
                if l.startswith('## Sample:'):
                    sname = l.strip().split()[-1]
                    print('Parsed sname', sname)
                if l.startswith('## Minimal completeness'):
                    m = re.match(r'## Minimal completeness: (?P<min_1x_pct>[\d.]+)% at 1x or (?P<min_5x_len>[\d.]+)bp at 5x', l)
                    assert m, l
                    min_1x_pct = m.group('min_1x_pct')
                    min_5x_len = m.group('min_5x_len')
                continue
            if l.startswith('#'):
                assert sname is not None, 'Could not parse sample name from the ## header'
                assert min_1x_pct is not None, 'Could not parse min_1x_pct from the ## header'
                assert min_5x_len is not None, 'Could not parse min_5x_len from the ## header'
                hdr_vals = l.strip().split("\t")  # #virus  size  depth  1x  5x  25x  significance
                continue
            values_dict = dict(zip(hdr_vals, l.strip().split("\t")))
            virus_name = values_dict['#virus']
            ave_depth = float(values_dict['depth'])
            pct_1x = float(values_dict['1x'])
            len_5x = int(float(values_dict['5x']) * int(values_dict['size']))
            significant = values_dict['significance'] == 'significant'
            viral_data.append((pct_1x, len_5x, ave_depth, virus_name, significant))
    viral_data.sort(reverse=True)
    there_some_hits = any(sign for p1x, l5x, d, v, sign in viral_data)
    if not there_some_hits:
        viral_data = [(p1x, l5x, d, v, sign) for p1x, l5x, d, v, sign in viral_data if sign]
    else:
        #TODO: use this for a separate oncoviral table:
        # showing all that significant, but at least 3 records even if nothing is significant
        viral_data = viral_data[:2]

    # To add into stats:
    metric = 'viral_content'
    data = {
        sname: {
            metric: "; ".join([v for i, (p1x, l5x, d, v, sign) in enumerate(viral_data)]) if there_some_hits else '-'
        }
    }
    header = {
        metric: dict(
            title='Viruses',
            description=(
                f'Detected viral sequences associated with cancer, per GDC database. '
                f'Format: x depth; 1x completeness (% of sequence covered at >1x), 5x length (base pairs covered at >5x). '
                f'Showing significant hits with 1x completeness at least {int(100 * pct_1x)}%, '
                f'or 5x length at least {len_5x}bp.'
            ),
        )
    }
    return data, header


def make_viral_integration_mqc_metrics(present_viruses, breakpoints_tsv):
    number_of_integration_sites_by_sample = defaultdict(int)
    with open(breakpoints_tsv) as f:
        for rec in csv.DictReader(f, delimiter='\t'):
            contig = rec['contig']
            sample_name = rec['sample']
            if contig not in present_viruses:  # human chromosome?
                number_of_integration_sites_by_sample[sample_name] += 1

    metric = 'viral_integration_sites'
    data = {sn: {metric: val} for sn, val in number_of_integration_sites_by_sample.items()}
    header = {
        metric: dict(
            title='Viral IS',
            description='The number of found viral integration sites',
            min=0,
            format='{:,.0f}',
        )
    }
    return data, header


def get_integration_sites_tsv_fn(wildcards):
    present_viruses = checkpoints.viral_content.get(**wildcards).output.present_viruses
    if open(present_viruses).read().strip():
        return present_viruses, rules.oncoviruses_breakpoints_tsv.output.tsv
    else:
        return present_viruses


rule oncoviral_multiqc:
    input:
        prioritized_tsv = '{batch}/oncoviruses/prioritized_oncoviruses.tsv',
        present_viruses = 'work/{batch}/oncoviruses/present_viruses.txt',
        wait_for_integration_sites = get_integration_sites_tsv_fn,
    output:
        data_yml = 'work/{batch}/oncoviruses/{batch}_oncoviruses_stats_data.yml',
        header_yml = 'work/{batch}/oncoviruses/{batch}_oncoviruses_stats_header.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name,
        breakpoints_tsv = 'work/{batch}/oncoviruses/oncoviral_breakpoints.tsv',
    run:
        data, header = make_oncoviral_mqc_metric(input.prioritized_tsv)
        headers = [header]

        with open(input.present_viruses) as f:
            present_viruses = [v for v in f.read().strip().split(',') if v]
        if present_viruses:
            integration_datas, integration_headers = \
                make_viral_integration_mqc_metrics(present_viruses, params.breakpoints_tsv)
            data = update_dict(data, integration_datas)
            headers.append(integration_headers)

        with open(output.data_yml, 'w') as out:
            yaml.dump(dict(
                data = data
            ), out, default_flow_style=False)

        # Using json instead of yaml because a long header wraps and breaks MultiQC
        with open(output.header_yml, 'w') as out:
            yaml.dump(dict(
                sp = dict(
                    oncoviruses = dict(
                        fn_re = basename(output.data_yml).replace(wildcards.batch, '.*'),
                    )
                ),
                custom_data = dict(
                    oncoviruses = dict(
                        plot_type = 'generalstats',
                        pconfig = headers
                    )
                ),
            ), out, default_flow_style=False)


rule oncoviruses_per_batch:
    input:
        get_integration_sites_tsv_fn,
        expand(rules.oncoviral_multiqc.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/oncoviruses_{batch}.done'))


rule oncoviruses:
    input:
        expand(rules.oncoviruses_per_batch.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/oncoviruses.done'))
