#!/usr/bin/env python
import glob
import os
import platform
import sys
from os import rename
from os.path import isfile, join, dirname, abspath, basename, exists
import click
import subprocess

from ngs_utils.call_process import run_simple
from ngs_utils.conda_utils import secondary_conda_env
from ngs_utils.file_utils import verify_file, safe_mkdir, splitext_plus, adjust_path
from ngs_utils import logger
from ngs_utils.logger import info, critical, err
from ngs_utils.snakemake_utils import run_snakemake
from ngs_utils.utils import set_locale; set_locale()
from umccrise import package_path
from reference_data import api as refdata


@click.command()
@click.option('-o', 'output_dir', type=click.Path(), help='Output directory (will be created if does not exist)')
@click.option('-T', 'tumor_bam', type=click.Path(exists=True), required=True, help='Tumor BAM')
@click.option('-N', 'normal_bam', type=click.Path(exists=True), required=False, help='Normal BAM')
@click.option('-nn', 'normal_name', help='Normal sample name')
@click.option('-tn', 'tumor_name', help='Tumor sample name')

# Ref data
@click.option('-g', 'genome', default='hg38')
@click.option('--genomes', '--genomes-dir', 'input_genomes_url', help='Path to umccrise genomes data. '
              'Can be s3 or gds. Alternative to --gridss-ref-dir')
@click.option('--ref-fa', 'ref_fa', help='Path to reference fasta (e.g. hg38)')
@click.option('--viruses-fa', 'viruses_fa', help='Path to viral sequences bwa index prefix '
              '(e.g. gdc-viral.fa from bcbio bundle)')
@click.option('--repeat-masker-bed', 'repeat_masker_bed', help='bedops rmsk2bed BED file for genome')

# Somatic filtering (GRIPSS)
@click.option('--breakend-pon', 'breakend_pon', help='Single breakend pon bed file, for somatic filtering (GRIPSS)')
@click.option('--bp-pon', 'bp_pon', help='Paired breakpoint pon bedpe file, for somatic filtering (GRIPSS)')
@click.option('--bp-hotspots', 'bp_hotspots', help='Paired breakpoint hotspot bedpe file '
              '(typically known fusion pairs) for somatic filtering (GRIPSS)')
@click.option('--min-tumor-af', 'min_tumor_af', type=click.FLOAT,
              help='Min tumor allelic frequency [0.005], for somatic filtering (GRIPSS)')

# Snakemake params:
@click.option('-t', '--threads', '-j', '--jobs', '--cores', 'requested_cores', type=click.INT,
              help='Maximum number of cores to use at single time (works both for local and cluster runs)')
@click.option('--unlock', 'unlock', is_flag=True)
@click.option('-n', '--dryrun', 'dryrun', is_flag=True,
              help='Propagated to snakemake. Prints rules and commands to be run without actually '
                   'executing them.')

# Gridss options
@click.option('--maxcoverage', 'maxcoverage', help='maximum coverage. Regions with coverage in excess of'
              ' this are ignored.')
@click.option('--chunksize-mil', 'chunksize_mil', help='In case if GRIDSS has attempted to open too many files '
              'at once and the OS file handle limit has been reached, it\'s useful to increase the chunk size '
              'from 10 million to 50 million using --chunksize-mil 50')
@click.option('--jvm_heap', 'jvm_heap', help='size of JVM heap for assembly and variant calling (-Xmx java parameter). '
              'Defaults to 25g to ensure GRIDSS runs on all cloud instances with approximate 32gb memory.'
              'However the maximum value is 32 to avoid Java\'s use of Compressed Oops which will effectively reduce '
              'the memory available to GRIDSS Recommendations: at least 4GB + 2GB per thread.')
@click.option('--externalaligner', 'externalaligner', type=click.Choice('bwa', 'minimap2'),
              help='Aligner to use: bwa or minimap2. BWA is more accurate, '
                   'but minimap2 does not require index and works much faster')

def main(output_dir=None, tumor_bam=None, normal_bam=None, normal_name=None, tumor_name=None,
         genome=None, input_genomes_url=None, ref_fa=None, viruses_fa=None, repeat_masker_bed=None,
         breakend_pon=None, bp_pon=None, bp_hotspots=None, min_tumor_af=None,
         requested_cores=None, unlock=False, dryrun=False,
         maxcoverage=None, chunksize_mil=None, jvm_heap=None, externalaligner=None):

    conf = {}

    output_dir = output_dir or 'gridss_results'
    output_dir = safe_mkdir(abspath(output_dir))
    log_dir = safe_mkdir(join(output_dir, 'log'))
    logger.init(log_fpath_=join(log_dir, 'gridss.log'), save_previous=True)
    if isfile(join(output_dir, 'work', 'all.done')):
        run_simple('rm ' + join(output_dir, 'work', 'all.done'))
    conf['output_dir'] = adjust_path(output_dir)

    tumor_name = tumor_name or splitext_plus(basename(tumor_bam))[0]\
        .replace('-ready', '').replace('-sorted', '')
    if normal_bam:
        normal_name = normal_name or splitext_plus(basename(normal_bam))[0]\
            .replace('-ready', '').replace('-sorted', '')
        conf['normal_bam'] = verify_file(normal_bam, 'Normal BAM, -N option'),
        conf['normal_name'] = normal_name
    conf['tumor_bam'] = verify_file(tumor_bam, 'Tumor BAM, -T option')
    conf['tumor_name'] = tumor_name

    try:
        machine_cores = len(os.sched_getaffinity(0))
    except:
        machine_cores = 1
    cores = min(machine_cores, 8)
    if requested_cores:
        cores = min(cores, requested_cores)
    conf['cores'] = cores

    if maxcoverage:
        conf['maxcoverage'] = maxcoverage
    if chunksize_mil:
        conf['chunksize_mil'] = chunksize_mil
    if jvm_heap:
        conf['jvm_heap'] = jvm_heap
    if externalaligner:
        conf['externalaligner'] = externalaligner

    conf['genome'] = genome
    try:
        from reference_data import api as refdata
    except:
        pass
    else:
        # check reference_data can find the genomes dir, and error out if not
        genomes_dir = refdata.find_genomes_dir(input_genomes_url)
        if genomes_dir:
            conf['genomes_dir'] = genomes_dir

    if ref_fa:
        if not externalaligner == 'minimap2' and not verify_file(ref_fa + '.bwt'):
            log.critical(f'Please, index {ref_fa} using'
                         f'    bwa index {ref_fa}')
        if not verify_file(ref_fa + '.fai'):
            log.critical(f'Please, index {ref_fa} using'
                         f'    samtools faidx {ref_fa}')
        conf['ref_fa'] = ref_fa
    if viruses_fa:
        if not externalaligner == 'minimap2' and not verify_file(viruses_fa + '.bwt'):
            log.critical(f'Please, index {viruses_fa} using: '
                         f'    bwa index {viruses_fa}')
        if not verify_file(viruses_fa + '.fai'):
            log.critical(f'Please, index {viruses_fa} using '
                         f'    samtools faidx {viruses_fa}')
        dict_file = viruses_fa.replace('.fa', '.dict')
        if not verify_file(dict_file):
            log.critical(f'Please, index {viruses_fa} using: '
                         f'   samtools dict {viruses_fa} -o {dict_file}')
        img_file = viruses_fa + '.img'
        if not verify_file(img_file):
            log.critical(f'Please, create an img file for {viruses_fa} using:\n'
                         f'   gatk BwaMemIndexImageCreator -I  {viruses_fa} -O {img_file}')

        conf['viruses_fa'] = verify_file(viruses_fa)
    if repeat_masker_bed:
        conf['repeat_masker_bed'] = repeat_masker_bed
    if breakend_pon:
        conf['breakend_pon'] = breakend_pon
    if bp_pon:
        conf['bp_pon'] = bp_pon
    if bp_hotspots:
        conf['bp_hotspots'] = bp_hotspots
    if min_tumor_af:
        conf['min_tumor_af'] = min_tumor_af

    py_path = sys.executable  # e.g. /miniconda/envs/umccrise_hmf/bin/python
    env_path = dirname(dirname(py_path))  # e.g. /miniconda/envs/umccrise_hmf
    found = glob.glob(join(env_path, 'share/gridss-*/gridss.jar'))
    if not found:
        hmf_env_path = secondary_conda_env('hmf', is_critical=False)
        if hmf_env_path:
            found = glob.glob(join(hmf_env_path, 'share/gridss-*/gridss.jar'))
            if not found:
                critical('Cannot find gridss JAR. Make sure you ran `conda install -c bioconda gridss`')
        conf['gridss_env'] = hmf_env_path
    conf['gridss_jar'] = found[0]


    run_snakemake(join(package_path(), 'gridss', 'Snakefile'), conf, cores=cores,
                  output_dir=output_dir, unlock=unlock, dryrun=dryrun)


if __name__ == '__main__':
    main()
