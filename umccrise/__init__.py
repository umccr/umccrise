import glob
import math
import os
import sys
from os.path import join, abspath, dirname, isfile, basename, splitext
from os.path import isfile, join, dirname, abspath, isdir
from datetime import datetime

from ngs_utils.Sample import BaseBatch, BaseProject, BaseSample
from ngs_utils.file_utils import splitext_plus, verify_file, verify_dir, adjust_path
from ngs_utils.bcbio import BcbioProject
from ngs_utils.dragen import DragenProject
from ngs_utils.utils import flatten
from ngs_utils.logger import critical, info, debug, warn, error
from ngs_utils import logger as log
from ngs_utils import bam_utils
from ngs_utils import vcf_utils
from ngs_utils.utils import set_locale; set_locale()
from ngs_utils.file_utils import verify_file
from reference_data import api as refdata


def package_path():
    return dirname(abspath(__file__))


def get_sig_rmd_file():
    """ Returns path to sig.Rmd file - R-markdown source for mutational signature analysys.
        The file must be located at the same directory as the Snakefile and the patient_analysis module.
    """
    return verify_file(join(package_path(), 'rmd_files', 'sig.Rmd'), is_critical=True)


class CustomSample(BaseSample):
    def __init__(self, **kwargs):
        BaseSample.__init__(self, **kwargs)  # name, dirpath, work_dir, bam, vcf, phenotype, normal_match
        self.qc_files = []

class CustomBatch(BaseBatch):
    # super class defines bam, somatic_vcf, germline_vcf, structural_vcf, qc_files, genome_build
    pass

class CustomProject(BaseProject):
    def __init__(self, input_dir=None, include_samples=None, exclude_samples=None,
                 silent=False, genome_build=None, **kwargs):
        BaseProject.__init__(self, input_dir=input_dir, **kwargs)
        self.include_samples = include_samples
        self.exclude_samples = exclude_samples
        self.genome_build = genome_build

    def get_or_create_batch(self, tumor_name, normal_name=None):
        if self.exclude_samples and tumor_name in self.exclude_samples:
            return None
        if self.include_samples and tumor_name not in self.include_samples:
            return None
        if normal_name:
            if self.exclude_samples and normal_name in self.exclude_samples:
                return None
            if self.include_samples and normal_name not in self.include_samples:
                return None

        if tumor_name in self.batch_by_name:
            b = self.batch_by_name[tumor_name]
        else:
            b = CustomBatch(tumor_name)
            b.tumor = CustomSample(name=tumor_name)
            if normal_name:
                b.normal = CustomSample(name=normal_name)
            self.batch_by_name[tumor_name] = b
            self.samples.extend([b.tumor, b.normal])
        return b

    def add_file(self, fpath):
        if fpath.endswith('.bam'):
            sample_name = bam_utils.sample_name_from_bam(fpath)
            b = self.get_or_create_batch(sample_name)
            if b:
                b.tumor.bam = fpath

        if fpath.endswith('.vcf.gz'):
            sample_names = vcf_utils.get_sample_names(fpath)
            if len(sample_names) == 1:
                b = self.get_or_create_batch(sample_names[0])
                if b:
                    b.somatic_vcf = fpath
            if len(sample_names) == 2:
                tumor_sn, normal_sn = sample_names
                b = self.get_or_create_batch(tumor_sn, normal_sn)
                if b:
                    b.somatic_vcf = fpath



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
        log.is_debug = True

    input_paths = smconfig.get('input_paths', [abspath(os.getcwd())])
    if isinstance(input_paths, str):
        input_paths = input_paths.split(',')

    custom_run = None
    found_bcbio_or_dragen_runs = []

    log.is_silent = silent  # to avoid redundant logging in cluster sub-executions of the Snakefile

    for input_path in input_paths:
        # custom file
        if isfile(input_path):
            if custom_run is None:
                custom_run = CustomProject(
                    input_dir=adjust_path(os.getcwd()),
                    include_samples=include_names,
                    exclude_samples=exclude_names)
            custom_run.add_file(input_path)

        # dragen
        elif isdir(input_path) and glob.glob(join(input_path, '*-replay.json')):
            run = DragenProject(input_path,
                               include_samples=include_names,
                               exclude_samples=exclude_names)
            found_bcbio_or_dragen_runs.append(run)
        # bcbio
        elif isdir(input_path):
            run = BcbioProject(input_path,
                               include_samples=include_names,
                               exclude_samples=exclude_names,
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
            run.batch_by_name = {b.name + '__' + b.tumor.name: b for b in batches}
            found_bcbio_or_dragen_runs.append(run)

        else:
            error(f'Cannot find file or dir {input_path}')

    if custom_run is None and len(found_bcbio_or_dragen_runs) == 1:
        # only one dragen or bcbio project - can return it directly
        combined_run = found_bcbio_or_dragen_runs[0]

    else:
        # multiple bcbio or dragen projects - combining them into a CustomProject
        combined_run = custom_run or CustomProject()
        for run in found_bcbio_or_dragen_runs:
            combined_run.batch_by_name.update(run.batch_by_name)
            combined_run.samples.extend(run.samples)
            if combined_run.project_name is None:
                combined_run.project_name = run.project_name
            if combined_run.genome_build is None:
                combined_run.genome_build = run.genome_build
            else:
                assert combined_run.genome_build == run.genome_build

    if combined_run.genome_build is None:
        combined_run.genome_build = 'hg38'
    if combined_run.project_name is None:
        combined_run.project_name = 'umccrise_' + datetime.now().strftime("%Y_%m_%d")

    # Reference files
    if smconfig.get('genomes_dir'):
        refdata.find_genomes_dir(smconfig.get('genomes_dir'))

    log.is_silent = False

    return combined_run, combined_run.batch_by_name


def prep_resources(num_batches, num_samples, ncpus_requested=None, is_cluster=False, is_silent=False):
    """ Determines the number of cpus used by a job and the total number of cpus
        available to snakemake scheduler.
        :returns ncpus_per_batch, ncpus_per_sample, ncpus_available, ncpus_per_node=None
    """
    # Checking presets for known HPC clusters, otherwise assuming a for single-machine AWS or local run
    # and just taking the number of available CPUs:

    if is_cluster:
        # we are not resticted to one machine, so can submit many jobs and let the scheduler figure out the queue
        ncpus_available = ncpus_requested or refdata.ncpus_on_node
        ncpus_per_node = refdata.ncpus_on_node
        if ncpus_per_node:
            info(f'Number of CPUs on a cluster node: {ncpus_per_node}')
        ncpus_local = 1
    else:
        try:
            ncpus_on_this_machine = len(os.sched_getaffinity(0))
        except:
            ncpus_on_this_machine = os.cpu_count()
        if ncpus_on_this_machine:
            info(f'Number of CPUs on this machine : {ncpus_on_this_machine}')
        # scheduling is on Snakemake, so restricting to the number of available cpus on this machine
        ncpus_available = min(ncpus_on_this_machine, ncpus_requested or math.inf)
        ncpus_per_node = None
        ncpus_local = ncpus_on_this_machine

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
        if ncpus_local:
            info(f'Local CPUs: {ncpus_local}')
        if ncpus_requested:
            info(f'Total CPUs requested by `umccrise --cores`: {ncpus_requested}')
        info(f'The pipeline can use {ncpus_available} CPUs total.')
        info(f'Batches found: {num_batches}, using {ncpus_per_batch} cpus per batch.')
        info(f'Samples found: {num_samples}, using {ncpus_per_sample} cpus per sample.')

    return ncpus_per_batch, ncpus_per_sample, ncpus_available, ncpus_per_node, ncpus_local


def prep_stages(include_stages=None, exclude_stages=None):
    default_enabled = {
        'conpair',
        'structural',
        'somatic', 'germline', 'maf',
        'purple',
        'mosdepth', 'goleft', 'cacao',
        'pcgr', 'cpsr',
        'oncoviruses',
        'rmd',
        'multiqc',
    }
    default_disabled = {
        'microbiome',
        'immuno',
    }
    debug(f'include_stages: {include_stages}')
    debug(f'exclude_stages: {exclude_stages}')
    def rename_input_stages(stages):
        fixed_stages = set()
        for s in stages:
            if s == 'cancer_report':
                fixed_stages |= {'rmd'}
            elif s == 'sv':
                fixed_stages |= {'structural'}
            elif s == 'purple':
                fixed_stages |= {'structural', 'purple'}
            elif s == 'coverage':
                fixed_stages |= {'mosdepth', 'goleft', 'cacao'}
            elif s == 'small_variants':
                fixed_stages |= {'somatic', 'germline'}
            elif s == 'cpsr':
                fixed_stages |= {'germline', 'cpsr'}
            elif s == 'pcgr':
                fixed_stages |= {'somatic', 'pcgr'}
            elif s not in default_enabled | default_disabled:
                warn(f'Stage {s} is not recognised. Available: {default_enabled | default_disabled}')
            elif s == 'default':
                fixed_stages |= default_disabled
            else:
                fixed_stages |= {s}
        return fixed_stages

    include_stages = rename_input_stages(include_stages)
    exclude_stages = rename_input_stages(exclude_stages)

    selected_stages = (include_stages or default_enabled) - exclude_stages
    debug(f'selected_stages: {selected_stages}')
    return selected_stages


def get_purple_metric(purple_file, metric='purity'):
    """ Reading the value from somatic sample from Purple output
    """
    with open(purple_file) as f:
        header, values = f.read().split('\n')[:2]
    # #Purity  NormFactor  Score   DiploidProportion  Ploidy  Gender  Status  PolyclonalProportion  MinPurity  MaxPurity  MinPloidy  MaxPloidy  MinDiploidProportion  MaxDiploidProportion  Version  SomaticDeviation
    # 0.7200   1.0400      0.3027  0.8413             1.8611  FEMALE  NORMAL  0.0000                0.6600     0.7700     1.8508     1.8765     0.8241                0.8558                2.17     0.0006
    data = dict(zip(header.strip('#').split('\t'), values.split('\t')))
    purity = float(data[metric])
    return purity


def get_purity(purple_file, phenotype='tumor'):
    """ Reading purity from somatic sample from Purple output
        Assuming purity 100% for normal
    """
    purity = 1.0
    if phenotype == 'tumor':
        purity = get_purple_metric(purple_file, 'purity')
    purity = min(purity, 1.0)
    return purity


def get_ploidy(purple_file):
    return get_purple_metric(purple_file, metric='ploidy')








