import glob
import math
import os
import sys
from os.path import join, abspath, dirname, isfile, basename, splitext
from os.path import isfile, join, dirname, abspath, isdir
from datetime import datetime

from ngs_utils.Sample import BaseBatch, BaseProject, BaseSample
from ngs_utils.file_utils import splitext_plus, verify_file, verify_dir, adjust_path
from ngs_utils.bcbio import BcbioProject, BcbioBatch
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
    def __init__(self, name, **kwargs):
        BaseBatch.__init__(self, name, **kwargs)  # bam, somatic_vcf, germline_vcf, structural_vcf, qc_files, genome_build
        self.qc_files = []
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


def _clean_sn(snames):
    return [v for v in flatten([sn.split('__') for sn in snames])]  # support "batch__sample" notation


def prep_inputs(smconfig, silent=False):
    ###############################
    #### Parsing a bcbio or Dragen project ####
    ###############################

    # Parsing a bcbio or Dragen project and including/excluding samples
    include_samples_map = dict()
    if smconfig.get('sample'):
        for sname in str(smconfig.get('sample')).split(','):
            if ':' in sname:
                ori_name, new_name = sname.split(':')
                include_samples_map[ori_name] = new_name
            else:
                include_samples_map[sname] = sname
    include_names = [ori_name for ori_name, new_name in include_samples_map.items()]

    exclude_names = smconfig.get('exclude')
    if exclude_names:
        exclude_names = str(exclude_names).split(',')

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
                    include_samples_map=include_samples_map,
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

    # Batch objects index by tumor sample names
    batches = [b for b in combined_run.batch_by_name.values()
               if not b.is_germline() and b.tumor and b.normal]
    assert batches
    batch_by_name = dict()
    for b in batches:
        new_name = include_samples_map.get(b.name) or \
                   include_samples_map.get(b.tumor.name) or \
                   b.name + '__' + b.tumor.name
        batch_by_name[new_name] = b

    return combined_run, batch_by_name


def prep_stages(include_stages=None, exclude_stages=None, run=None):
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
    # if not all(b.germline_vcf for b in run.batch_by_name.values()):
    #     default_disabled |= {'germline'}
    # if not all(b.somatic_vcf for b in run.batch_by_name.values()):
    #     default_disabled |= {'somatic'}
    # if not all(b.sv_vcf for b in run.batch_by_name.values()):
    #     default_disabled |= {'structural'}

    debug(f'include_stages: {include_stages}')
    debug(f'exclude_stages: {exclude_stages}')
    def _rename_input_stages(stages):
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
            elif s == 'oncoviral':
                fixed_stages |= {'oncoviruses'}
            elif s not in default_enabled | default_disabled:
                critical(f'Stage "{s}" is not recognised. Available: {default_enabled | default_disabled}')
            elif s == 'default':
                fixed_stages |= default_disabled
            else:
                fixed_stages |= {s}
        return fixed_stages

    include_stages = _rename_input_stages(include_stages)
    exclude_stages = _rename_input_stages(exclude_stages)

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








