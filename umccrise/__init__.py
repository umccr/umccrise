import glob
import math
import os
import sys
import multiprocessing
from os.path import join, abspath, dirname, isfile, basename, splitext
from ngs_utils.file_utils import splitext_plus, verify_file, verify_dir, adjust_path
from ngs_utils.bcbio import BcbioProject
from ngs_utils.dragen import DragenProject
from ngs_utils.utils import flatten
from ngs_utils.logger import critical, info, debug, warn
from ngs_utils import logger as ngs_utils_logger
from hpc_utils import hpc
from ngs_utils.utils import set_locale; set_locale()
from os.path import isfile, join, dirname, abspath
from ngs_utils.file_utils import verify_file


def package_path():
    return dirname(abspath(__file__))


def get_sig_rmd_file():
    """ Returns path to sig.Rmd file - R-markdown source for mutational signature analysys.
        The file must be located at the same directory as the Snakefile and the patient_analysis module.
    """
    return verify_file(join(package_path(), 'rmd_files', 'sig.Rmd'), is_critical=True)


def prep_inputs(smconfig, silent=False):
    ###############################
    #### Parsing a bcbio or Dragen project ####
    ###############################

    # Parsing a bcbio or Dragen project and including/excluding samples
    include_names = smconfig.get('batch') or smconfig.get('sample')
    if include_names:
        include_names = str(include_names).split(',')
        include_names = [v for v in flatten([sn.split('__') for sn in include_names])]  # support "batch__sample" notation
    exclude_names = smconfig.get('exclude')
    if exclude_names:
        exclude_names = str(exclude_names).split(',')
        exclude_names = [v for v in flatten([sn.split('__') for sn in exclude_names])]  # support "batch__sample" notation

    if smconfig.get('debug', 'no') == 'yes':
        ngs_utils_logger.is_debug = True

    ngs_utils_logger.is_silent = silent  # to avoid redundant logging in cluster sub-executions of the Snakefile
    input_path = smconfig.get('input_project', abspath(os.getcwd()))
    if glob.glob(join(input_path, '*-replay.json')):
        run = DragenProject(input_path)
        batch_by_name = run.batch_by_name

    else:
        run = BcbioProject(input_path,
                           exclude_samples=exclude_names,
                           include_samples=include_names,
                           silent=True)
        run.project_name = splitext(basename(run.bcbio_yaml_fpath))[0]

        if len(run.batch_by_name) == 0:
            if exclude_names:
                critical(f'Error: no samples left with the exclusion of batch/sample name(s): {", ".join(exclude_names)}.'
                         f'Check yaml file for available options: {run.bcbio_yaml_fpath}.')
            if include_names:
                critical(f'Error: could not find a batch or a sample with the name(s): {", ".join(include_names)}. '
                         f'Check yaml file for available options: {run.bcbio_yaml_fpath}')
            critical(f'Error: could not parse any batch or samples in the bcbio project. '
                     f'Please check the bcbio yaml file: {run.bcbio_yaml_fpath}')

        # Batch objects index by tumor sample names
        batches = [b for b in run.batch_by_name.values() if not b.is_germline() and b.tumor and b.normal]
        assert batches
        batch_by_name = {b.name + '__' + b.tumor.name: b for b in batches}

    # Reference files
    if smconfig.get('genomes_dir'):
        hpc.set_genomes_dir(smconfig.get('genomes_dir'))

    ngs_utils_logger.is_silent = False
    return run, batch_by_name


def prep_resources(num_batches, num_samples, ncpus_requested=None, is_cluster=False, is_silent=False):
    """ Determines the number of cpus used by a job and the total number of cpus
        available to snakemake scheduler.
        :returns ncpus_per_batch, ncpus_per_sample, ncpus_available, ncpus_per_node=None
    """
    # Checking presets for known HPC clusters, otherwise assuming a for single-machine AWS or local run
    # and just taking the number of available CPUs:
    ncpus_on_a_machine = hpc.ncpus_on_node or multiprocessing.cpu_count() or 1
    if is_cluster:
        # we are not resticted to one machine, so can submit many jobs and let the scheduler figure out the queue
        ncpus_available = ncpus_requested or 32
        ncpus_per_node = ncpus_on_a_machine
    else:
        # scheduling is on Snakemake, so restricting to the number of available cpus on a machine
        ncpus_available = min(ncpus_on_a_machine, ncpus_requested or math.inf)
        ncpus_per_node = None

    ncpus_per_batch = max(1, ncpus_available // num_batches)
    ncpus_per_sample = max(1, ncpus_available // num_samples)

    def adjust_ncpus_per_job(ncpus, max_ncpus_per_job=10, msg=''):
        """ Adjusting the number of cpus to a number below <max_ncpus_per_job>.
            Say, if we have more than 20 cpus on a node and only 1 batch, we should adjust
            to use only half of that for a batch, so that 2 different jobs (say, AMBER and COBALT)
            can be run in parallel, because using 20 cpus per one job is a waste.
        """
        if ncpus > max_ncpus_per_job:
            # new_ncpus = ncpus
            factor = math.ceil(ncpus / max_ncpus_per_job)
            new_ncpus = ncpus // factor
            # while True:
            #     factor += 1
            #     new_ncpus = ncpus // factor
            #     print(f'ncpus: {ncpus}, factor: {factor}, new_ncpus: {new_ncpus}')
            #     if new_ncpus < max_ncpus_per_job:
            #         print(f'breaking')
            #         break
            if not is_silent:
                info((msg if msg else 'The number of cpus per batch is ') + f'{ncpus} >{max_ncpus_per_job}. '
                     f'This is usually wasteful, so we are adjusting it '
                     f'to the number <={max_ncpus_per_job}: {new_ncpus} = {ncpus} // {factor}, so '
                     f'{factor} different rules can be run in parallel (say, AMBER and COBALT '
                     f'at the same time).')
            ncpus = new_ncpus
        return ncpus

    ncpus_per_batch = adjust_ncpus_per_job(ncpus_per_batch, max_ncpus_per_job=14, msg=
        f'The number of cpus available is {ncpus_available}, and the number of batches is {num_batches}, '
        f'so the number of cpus per batch would be ')
    ncpus_per_sample = adjust_ncpus_per_job(ncpus_per_sample, max_ncpus_per_job=14, msg=
        f'The number of cpus available is {ncpus_available}, and the number of samples is {num_samples}, '
        f'so the number of cpus per sample would be ')

    if not is_silent:
        info(f'Final number of CPUs per machine: {ncpus_on_a_machine}')
        if ncpus_requested:
            info(f'Total CPUs requested by `umccrise -t`: {ncpus_requested}')
        info(f'The pipeline can use {ncpus_available} CPUs total.')
        info(f'Batches found: {num_batches}, using {ncpus_per_batch} cpus per batch.')
        info(f'Samples found: {num_samples}, using {ncpus_per_sample} cpus per sample.')

    return ncpus_per_batch, ncpus_per_sample, ncpus_available, ncpus_per_node
