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

# Snakemake params:
@click.option('-t', '--threads', '-j', '--jobs', '--cores', 'requested_cores', type=click.INT,
              help='Maximum number of cores to use at single time (works both for local and cluster runs)')
@click.option('--unlock', 'unlock', is_flag=True)
@click.option('-n', '--dryrun', 'dryrun', is_flag=True,
              help='Propagated to snakemake. Prints rules and commands to be run without actually '
                   'executing them.')

def main(output_dir=None, tumor_bam=None, normal_bam=None, normal_name=None, tumor_name=None,
         genome=None, input_genomes_url=None, ref_fa=None, viruses_fa=None,
         requested_cores=None, unlock=False, dryrun=False):

    conf = {}

    output_dir = output_dir or 'gridss_results'
    output_dir = safe_mkdir(abspath(output_dir))
    logger.init(log_fpath_=join(output_dir, 'gridss.log'), save_previous=True)
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
        cores = len(os.sched_getaffinity(0))
    except:
        cores = 1
    if requested_cores:
        cores = min(requested_cores, cores)
    conf['cores'] = cores

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
        if not verify_file(ref_fa + '.bwt'):
            log.critical(f'Please, index {ref_fa} with `bwa index {ref_fa}`')
        if not verify_file(ref_fa + '.fai'):
            log.critical(f'Please, index {ref_fa} with `samtools faidx {ref_fa}`')
        conf['ref_fa'] = ref_fa
    if viruses_fa:
        if not verify_file(viruses_fa + '.bwt'):
            log.critical(f'Please, index {viruses_fa} with `bwa index {viruses_fa}`')
        if not verify_file(viruses_fa + '.fai'):
            log.critical(f'Please, index {viruses_fa} with `samtools faidx {viruses_fa}`')
        dict_file = viruses_fa.replace('.fa', '.dict')
        if not verify_file(dict_file):
            log.critical(f'Please, index {viruses_fa} with `samtools dict {viruses_fa} -o {dict_file}`')
        conf['viruses_fa'] = verify_file(viruses_fa)

    py_path = sys.executable  # e.g. /miniconda/envs/umccrise_hmf/bin/python
    env_path = dirname(dirname(py_path))  # e.g. /miniconda/envs/umccrise_hmf
    found = glob.glob(join(env_path, 'share/gridss-*/gridss.jar'))
    if not found:
        hmf_env_path = secondary_conda_env('hmf', is_critical=False)
        if hmf_env_path:
            found = glob.glob(join(hmf_env_path, 'share/gridss-*/gridss.jar'))
            if not found:
                critical('Cannot find gridss JAR. Make sure you ran `conda install -c bioconda gridss`')
    conf['gridss_jar'] = found[0],

    run_snakemake(join(package_path(), 'gridss', 'Snakefile'), conf, cores=cores,
                  output_dir=output_dir, unlock=unlock, dryrun=dryrun)


if __name__ == '__main__':
    main()
