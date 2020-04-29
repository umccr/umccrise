import json
from os.path import join, basename
import yaml
from hpc_utils import hpc
from ngs_utils.file_utils import verify_file
from cyvcf2 import VCF


localrules: oncoviruses, oncoviruses_per_batch


checkpoint viral_content:
    input:
        tumor_bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
    output:
        prioritized_tsv = '{batch}/oncoviruses/prioritized_oncoviruses.tsv',
        present_viruses = 'work/{batch}/oncoviruses/present_viruses.txt',
    params:
        genomes_dir = hpc.genomes_dir,
        work_dir = 'work/{batch}/oncoviruses',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
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
        genomes_dir = hpc.genomes_dir,
        work_dir = 'work/{batch}/oncoviruses',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
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
              '--genomes-dir {params.genomes_dir} {params.unlock_opt} -v $(cat {input.significant_viruses})')
        shell('cp {params.work_dir}/breakpoints.vcf.gz {output.breakpoints_vcf}')


# TODO: count viruses, breakpoints and genes for MultiQC report
def make_oncoviral_mqc_metric(prioritized_oncoviruses_tsv):
    viral_data = []
    with open(prioritized_oncoviruses_tsv) as f:
        ## Viral sequences (from {ONCOVIRAL_SOURCE_URL}) found in unmapped reads
        ## Sample: {SAMPLE}
        ## Significant completeness: {MIN_SIGNIFICANT_COMPLETENESS}
        ## Significant coverage: {COMPLETENESS_THRESHOLD}
        #virus\tsize\tdepth\t1x\t5x\t25x\tsignificance
        for l in f:
            if l.startswith('##'):
                if l.startswith('## Sample'):
                    sample_name = l.strip().split()[-1]
                if l.startswith('## Significant completeness'):
                    significant_completeness = float(l.strip().split()[-1])
                if l.startswith('## Significant coverage'):
                    significant_coverage = l.strip().split()[-1]
                continue
            if l.startswith('#'):
                headers = l.strip().split("\t")  # #virus  size    depth   1x      5x      25x
                continue
            values_dict = dict(zip(headers, l.strip().split("\t")))
            virus_name = values_dict['#virus']
            ave_depth = float(values_dict['depth'])
            try:
                completeness = float(values_dict[significant_coverage])
            except:
                print(significant_coverage, values_dict)
                raise
            viral_data.append((completeness, ave_depth, virus_name))
    viral_data.sort(reverse=True)
    there_some_hits = any(c >= significant_completeness for c, d, v in viral_data)
    if not there_some_hits:
        # showing all that significant, but at least 3 records even if nothing is significant
        viral_data = viral_data[:2]
    else:
        viral_data = [(c, d, v) for c, d, v in viral_data if c >= significant_completeness]

    # To add into stats:
    data = {
        sample_name: {
            "viral_content": "; ".join([v for i, (c, d, v) in enumerate(viral_data)]) if there_some_hits else '-'
        }
    }
    header = {
        'title': 'Viral',
        'description': (
            f'Detected viral sequences associated with cancer, per GDC database. '
            f'Format: x depth; completeness (% of sequence covered at >{significant_coverage}). '
            f'Showing significant hits with at least {int(100 * significant_completeness)}% completeness.'
        ),
    }
    return data, header


rule oncoviral_multiqc:
    input:
        prioritized_tsv = '{batch}/oncoviruses/prioritized_oncoviruses.tsv'
    output:
        data_yml = 'work/{batch}/oncoviruses/{batch}_oncoviruses_stats_data.yml',
        header_yml = 'work/{batch}/oncoviruses/{batch}_oncoviruses_stats_header.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name
    group: "oncoviruses"
    run:
        data, header = make_oncoviral_mqc_metric(input.prioritized_tsv)

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
                        pconfig = [dict(
                            viral_content = header,
                        )]
                    )
                ),
            ), out, default_flow_style=False)


def parse_info_field(rec, name):
    val = rec.INFO.get(name)
    if val is None:
        return ''
    elif isinstance(val, float) or isinstance(val, int) or isinstance(val, bool) or isinstance(val, str):
        return str(val)
    else:
        return ','.join(map(str, val))


# Produce a TSV file for further analysis in Rmd
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
    group: "oncoviruses"
    run:
        sample_name = batch_by_name[wildcards.batch].tumor.name
        viruses = open(input.present_viruses).read().split(',')
        with open(output.tsv, 'w') as out:
            header = ['sample', 'contig', 'start', 'end', 'svtype',
                      'PAIR_COUNT', 'Genes',  'ID', 'MATEID',]
            out.write('\t'.join(header) + '\n')
            for rec in VCF(input.vcf):
                viral_genes = rec.INFO.get('ViralGenes')
                host_genes = rec.INFO.get('GenesWithin100kb')
                data = [sample_name, rec.CHROM, rec.POS, rec.INFO.get('END', ''),
                        rec.INFO['SVTYPE'],
                        parse_info_field(rec, 'PAIR_COUNT'),
                        viral_genes if rec.CHROM in viruses else host_genes,
                        rec.ID,
                        parse_info_field(rec, 'MATEID'),
                        ]
                out.write('\t'.join(map(str, data)) + '\n')


def get_integration_sites_tsv_fn(wildcards):
    present_viruses = checkpoints.viral_content.get(**wildcards).output.present_viruses
    if open(present_viruses).read().strip():
        return present_viruses, rules.oncoviruses_breakpoints_tsv.output.tsv
    else:
        return present_viruses


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
