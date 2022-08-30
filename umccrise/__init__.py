import subprocess
import glob
import math
import os
import pathlib
import re
import sys
from os.path import join, abspath, dirname, isfile, basename, splitext
from os.path import isfile, join, dirname, abspath, isdir
from datetime import datetime
import csv

import pysam

from ngs_utils.Sample import BaseBatch, BaseProject, BaseSample
from ngs_utils.file_utils import splitext_plus, verify_file, verify_dir, adjust_path
from ngs_utils.bcbio import BcbioProject, BcbioBatch, detect_bcbio_dir
from ngs_utils.dragen import DragenProject
from ngs_utils.utils import flatten
from ngs_utils.logger import critical, info, debug, warn, error
from ngs_utils import logger as log
from ngs_utils import bam_utils
from ngs_utils import vcf_utils
from ngs_utils.utils import set_locale; set_locale()
from ngs_utils.file_utils import verify_file, verify_obj_by_path
from reference_data import api as refdata
from ngs_utils.reference_data import get_all_genes_bed


def package_path():
    return dirname(abspath(__file__))


def get_sig_rmd_file():
    """ Returns path to sig.Rmd file - R-markdown source for mutational signature analysys.
        The file must be located at the same directory as the Snakefile and the patient_analysis module.
    """
    return verify_file(join(package_path(), 'rmd_files', 'sig.Rmd'), is_critical=True)


class CustomSample(BaseSample):
    def __init__(self, assay=None, **kwargs):
        BaseSample.__init__(self, **kwargs)  # name, dirpath, work_dir, bam, vcf, phenotype, normal_match
        self.qc_files = []
        self.assay = assay

class CustomBatch(BaseBatch):
    def __init__(self, name, **kwargs):
        BaseBatch.__init__(self, name, **kwargs)  # bam, somatic_vcf, germline_vcf, structural_vcf, qc_files, genome_build
        self.qc_files = []
    pass

class CustomProject(BaseProject):
    def __init__(self, input_dir=None,
                 include_samples=None, exclude_samples=None,
                 silent=False, genome_build=None, **kwargs):
        BaseProject.__init__(self, input_dir=input_dir, **kwargs)
        self.include_samples = include_samples
        self.exclude_samples = exclude_samples
        self.genome_build = genome_build
        self.input_tsv_fpaths = []
        self._parsed_bcbio_projects_by_path = dict()

    def _load_bcbio_project(self, bcbio_project_path):
        proj = self._parsed_bcbio_projects_by_path.get(bcbio_project_path)
        if not proj:
            info(f'Loading project {bcbio_project_path}')
            proj = BcbioProject(bcbio_project_path, silent=True)
            self._parsed_bcbio_projects_by_path[bcbio_project_path] = proj
        return proj

    def add_batch(self, entry, base_path):
        if self.exclude_samples and entry['sample'] in self.exclude_samples:
            return None
        if self.include_samples and entry['sample'] not in self.include_samples:
            return None

        b = CustomBatch(name=entry['sample'],
                        somatic_vcf=entry.get('somatic_vcf'),
                        germline_vcf=entry.get('germline_vcf'),
                        sv_vcf=entry.get('sv_vcf'))
        self.batch_by_name[entry['sample']] = b

        def _full_path(path):
            if not path or path == '.':
                return None
            if path.startswith('/'):
                return verify_obj_by_path(path, is_critical=True)
            else:
                path = join(base_path, path)
                return verify_obj_by_path(path, is_critical=True)

        # Adding RNA project for neoantigens
        rna_bcbio_path = _full_path(entry.get('rna_bcbio'))
        rna_sname = entry.get('rna_sample')
        rna_sample = None
        if rna_bcbio_path:
            if not rna_sname:
                critical(f'rna_sample must be provided along with rna_bcbio '
                         f'(for sample {entry["sample"]})')
            rna_bcbio_project = self._load_bcbio_project(rna_bcbio_path)
            rna_samples = [s for s in rna_bcbio_project.samples
                          if s.name == rna_sname]
            if not rna_samples:
                critical(f'Could not find RNA sample {rna_sname} in {rna_bcbio_path}')
            rna_sample = rna_samples[0]

        # Adding BAMs
        rna_bam = _full_path(entry.get('rna'))
        if not rna_bam and rna_sample:
            rna_bam = rna_sample.bam
        wgs_bam = _full_path(entry.get('wgs'))
        normal_bam = _full_path(entry.get('normal'))
        exome_bam = _full_path(entry.get('exome'))
        exome_normal_bam = _full_path(entry.get('exome_normal'))

        if wgs_bam:
            wgs_s = CustomSample(
                name=entry['sample'],
                rgid=bam_utils.sample_name_from_bam(wgs_bam),
                phenotype='tumor',
                assay='wgs',
                bam=wgs_bam)
            self.samples.append(wgs_s)
            b.add_tumor(wgs_s)

            if normal_bam:
                wgs_normal_s = CustomSample(
                    name=entry['sample'] + '_normal',
                    rgid=bam_utils.sample_name_from_bam(normal_bam),
                    phenotype='normal',
                    assay='wgs',
                    bam=normal_bam)
                self.samples.append(wgs_normal_s)
                wgs_s.normal_matches.append(wgs_normal_s)
                b.add_normal(wgs_normal_s)

        if exome_bam:
            exome_s = CustomSample(
                name=entry['sample'] + '_exome',
                rgid=bam_utils.sample_name_from_bam(exome_bam),
                phenotype='tumor',
                assay='exome',
                bam=exome_bam)
            self.samples.append(exome_s)
            b.add_tumor(exome_s)

            if exome_normal_bam:
                exome_normal_s = CustomSample(
                    name=entry['sample'] + '_normal_exome',
                    rgid=bam_utils.sample_name_from_bam(exome_normal_bam),
                    phenotype='normal',
                    assay='exome',
                    bam=exome_normal_bam)
                self.samples.append(exome_normal_s)
                exome_s.normal_matches.append(exome_normal_s)
                b.add_normal(exome_normal_s)

        if rna_sample:
            self.samples.append(rna_sample)
            b.add_rna_sample(rna_sample)
        elif rna_bam:
            rna_s = CustomSample(
                name=entry['sample'] + '_rna',
                rgid=bam_utils.sample_name_from_bam(rna_bam),
                phenotype='normal',
                assay='rna',
                bam=rna_bam)
            self.samples.append(rna_s)
            b.add_rna_sample(rna_s)

        return b

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
            b.tumors.append(CustomSample(name=tumor_name))
            if normal_name:
                b.normals.append(CustomSample(name=normal_name))
            self.batch_by_name[tumor_name] = b
            self.samples.extend(b.tumors + b.normals)
        return b

    def add_file(self, fpath):
        if fpath.endswith('.bam'):
            sample_name = bam_utils.sample_name_from_bam(fpath)
            b = self.get_or_create_batch(sample_name)
            if b:
                b.tumors[0].bam = fpath

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

    # Handle inputs provided by positional arguments
    input_paths_raw = smconfig.get('input_paths')
    if input_paths_raw:
        if isinstance(input_paths_raw, list):
            input_paths_strs = input_paths_raw
        elif isinstance(input_paths_raw, str):
            input_paths_strs = input_paths_raw.split(',')
        else:
            assert False
        input_paths = [pathlib.Path(p) for p in input_paths_strs]
    else:
        input_paths = list()

    # Handle inputs provided by named arguments. DRAGEN directories are handled separately below.
    dragen_directories = list()
    bcbio_directories = smconfig.get('bcbio_dirs', list())
    custom_tsvs = smconfig.get('custom_tsv_fps', list())

    found_bcbio_or_dragen_runs = []
    custom_run = CustomProject(
        input_dir=adjust_path(os.getcwd()),
        include_samples=include_names,
        exclude_samples=exclude_names
    )

    log.is_silent = silent  # to avoid redundant logging in cluster sub-executions of the Snakefile

    # Recursively search for bcbio and DRAGEN output directories. We do *not* recurse into
    # bcbio/DRAGEN output directories i.e. nested outputs are ignored
    #
    # A recursive search is used to simplify and generalise DRAGEN output discovery by not relying
    # on file name patterns, which are subject to change. See DRAGEN output structuring (20211108):
    #   * DRAGEN tumor/normal:
    #       - <root_output_dir>/<subject_id>_dragen_somatic/<output_files>
    #   * DRAGEN normal:
    #       - <root_output_dir>/<subject_id>_dragen_germline/<output_files>
    # The DRAGEN tumor/normal and normal output directories may be staged under a single root
    # directory.
    #
    # Input custom TSVs must be explicitly provided on the command line.
    #
    # Currently to generate DRAGEN tumor/normal data, two DRAGEN runs are done: (1) tumor/normal,
    # and (2) normal. Hence there the two DRAGEN output directories must be paired prior to
    # launching the main umccrise Snakemake workflow.
    #
    # Initialisation of SampleProject deriviative instances is done below to allow pairing of
    # DRAGEN output directories.

    # First take all explicitly set *files* from the input. These will be custom TSVs, BAMs, and/or
    # VCFs.
    provided_files = list()
    directories = list()
    invalid_inputs = list()
    for path in input_paths:
        if path.name.endswith('.tsv'):
            custom_tsvs.append(path)
        elif any(path.name.endswith(fext) for fext in ('.bam', '.vcf.gz')):
            provided_files.append(path)
        elif path.is_dir():
            directories.append(path)
        else:
            invalid_inputs.append(path)
    if invalid_inputs:
        critical(f'received bad command line inputs: {", ".join(invalid_inputs)}')

    # Discover bcbio and DRAGEN directories
    for directory in directories:
        discovered_directories = recursively_search_input(directory)
        if not any(dl for dl in discovered_directories.values()):
            critical(f'did not find any valid inputs in {directory}')
        bcbio_directories.extend(discovered_directories['bcbio'])
        dragen_directories.extend(discovered_directories['dragen'])

    # Pair DRAGEN directories if given by input_paths otherwise manually create and collect
    # required information to initialise a DragenProject from the provided paired DRAGEN input
    # directories (i.e. from --dragen_somatic_dir and --dragen_germline_dir)
    if any(arg in smconfig for arg in ('dragen_somatic_dir', 'dragen_germline_dir')):
        assert len(dragen_directories) == 0
        dragen_run_sets = create_dragen_paired_directories_from_config(smconfig)
    else:
        dragen_run_sets = pair_dragen_directories(dragen_directories)

    # Create input project classes
    # Custom TSVs
    for custom_tsv in custom_tsvs:
        custom_run.input_tsv_fpaths.append(custom_tsv)
        with open(custom_tsv) as f:
            for entry in csv.DictReader(f, delimiter='\t'):
                b = custom_run.add_batch(entry, dirname(custom_tsv))
                if b:
                    include_samples_map[b.name] = b.name
    # BAMs/VCFs provided directly on command line
    for provided_file in provided_files:
        custom_run.add_file(provided_file)
    # Paired DRAGEN directories
    for dragen_run_set in dragen_run_sets:
        # Include or exclude run based on available identifiers
        identifiers = [
            dragen_run_set['subject_id'],
            dragen_run_set['normal_run']['normal'],
            dragen_run_set['tumor_normal_run']['tumor'],
            dragen_run_set['tumor_normal_run']['normal'],
        ]
        if include_names and all(ident not in include_names for ident in identifiers):
            continue
        if exclude_names and any(ident in exclude_names for ident in identifiers):
            continue
        # Initialise project
        run = DragenProject(
            dragen_run_set,
            include_samples=include_names,
            exclude_samples=exclude_names
        )
        found_bcbio_or_dragen_runs.append(run)
    # bcbio directories
    for bcbio_directory in bcbio_directories:
        run = BcbioProject(
            bcbio_directory,
            include_samples=include_names,
            exclude_samples=exclude_names,
            silent=True
        )
        run.project_name = splitext(basename(run.bcbio_yaml_fpath))[0]
        if not run.batch_by_name:
            continue
        found_bcbio_or_dragen_runs.append(run)

    # NOTE(SW): previously if multiple Batches were supplied, a CustomProject was initialised to
    # combine them into a single project. However, the CustomProject class is not compatible with
    # the MultiQC stage; the process of collecting QC files is input-type specific and implemented
    # in corresponding Project/Batch classes (e.g. DragenProject, BcbioBatch) but the
    # CustomProject/CustomBatch only implements a generalised QC files collection function.
    #
    # Consequently, CustomProject and DragenProject cause issues with the MultiQC report, and so
    # multiple Batches have been disallowed so that well-formed MultiQC reports can be generated.
    # This limitation will likely be lifted in the future.
    #
    # Currently umccrise accepts either one DRAGEN run, one bcbio batch, or one 'custom' run
    if len(custom_run.samples) == 0 and len(found_bcbio_or_dragen_runs) == 1:
        run = found_bcbio_or_dragen_runs[0]
    elif len(custom_run.samples) == 1 and len(found_bcbio_or_dragen_runs) == 0:
        run = custom_run
    else:
        runs_count = len([*custom_run.samples, *found_bcbio_or_dragen_runs])
        if runs_count == 0:
            critical('no batches found')
        else:
            critical(
                f'discovered {runs_count} batches but umccrise currently only accepts '
                'a single batch'
            )

    # NOTE(SW): we are not currently supporting backwards compatibility for bcbio data and hence
    # reject bcbio input. Known required changes to achieve bcbio backwards compatibility are:
    #   - updating MultiQC mosdepth module to display sample names consistently with other modules
    #   - regenerating QC reference bcbio data to propogate changes
    if bcbio_directories:
        critical('the current version of umccrise does not support bcbio input')

    if run.genome_build is None:
        run.genome_build = 'hg38'
    if run.project_name is None:
        run.project_name = 'umccrise_' + datetime.now().strftime("%Y_%m_%d")

    # Reference files
    if smconfig.get('genomes_dir'):
        refdata.find_genomes_dir(smconfig.get('genomes_dir'))

    log.is_silent = False

    # Ensure at least one batch has been discovered
    batches = [b for b in run.batch_by_name.values()
               if not b.is_germline() and b.tumors
               and b.normals]
    if not batches:
        msgs = list()
        msgs.append('No batches discovered')
        if include_names:
            msgs.append(f'Include batch/sample name list: {", ".join(include_names)}')
        if exclude_names:
            msgs.append(f'Exclude batch/sample name list: {", ".join(exclude_names)}')
        critical('. '.join(msgs))
    # Batch objects index by tumor sample names
    batch_by_name = dict()
    for b in batches:
        new_name = include_samples_map.get(b.name) or \
                   include_samples_map.get(b.tumors[0].name) or \
                   b.name + '__' + b.tumors[0].name
        batch_by_name[new_name] = b

    return run, batch_by_name


def prep_stages(run, include_stages=None, exclude_stages=None):
    default_enabled = {
        'conpair',
        'structural',
        'somatic', 'germline', 'maf',
        'purple',
        'mosdepth', 'goleft', 'cacao',
        'pcgr', 'cpsr',
        'oncoviruses',
        'cancer_report',
        'multiqc',
        'combined_multiqc',
    }
    default_disabled = {
        'microbiome',
        'neoantigens',
        'peddy',
        'haplotype_caller',
    }
    # if not all(b.germline_vcf for b in run.batch_by_name.values()):
    #     default_disabled |= {'germline'}
    # if not all(b.somatic_vcf for b in run.batch_by_name.values()):
    #     default_disabled |= {'somatic'}
    # if not all(b.sv_vcf for b in run.batch_by_name.values()):
    #     default_disabled |= {'structural'}

    # For non-bcbio runs we want to run BCFtools and SAMtools stats on input BAMs by default
    if not isinstance(run, BcbioProject):
        default_enabled |= {'samtools_stats', 'bcftools_stats'}

    debug(f'include_stages: {include_stages}')
    debug(f'exclude_stages: {exclude_stages}')
    def _rename_input_stages(stages):
        fixed_stages = set()
        for s in stages:
            if s in {'rmd', 'cancer_report'}:
                fixed_stages |= {'cancer_report', 'oncoviruses'}
            elif s in {'sv', 'structural'}:
                fixed_stages |= {'structural'}
            elif s in {'purple', 'cnv'}:
                fixed_stages |= {'structural', 'purple'}
            elif s in {'coverage'}:
                fixed_stages |= {'mosdepth', 'goleft', 'cacao'}
                if not isinstance(run, BcbioProject):
                    fixed_stages |= {'samtools_stats'}
            elif s == 'small_variants':
                fixed_stages |= {'somatic', 'germline'}
            elif s == 'cpsr':
                fixed_stages |= {'germline', 'cpsr'}
            elif s == 'pcgr':
                fixed_stages |= {'somatic', 'pcgr'}
            elif s in {'oncoviral', 'oviraptor', 'oncoviruses', 'viruses', 'viral'}:
                fixed_stages |= {'oncoviruses'}
            elif s in {'nag', 'immuno', 'neoantigens'}:
                fixed_stages |= {'neoantigens'}
            elif s == 'multiqc':
                fixed_stages |= {'multiqc', 'purple', 'conpair', 'somatic', 'germline',
                                 'oncoviruses', 'mosdepth'}
                if not isinstance(run, BcbioProject):
                    fixed_stages |= {'samtools_stats', 'bcftools_stats'}
            elif s == 'default':
                fixed_stages |= default_enabled
            elif s not in default_enabled | default_disabled:
                critical(f'Stage "{s}" is not recognised. Available: {default_enabled | default_disabled}')
            else:
                fixed_stages |= {s}
        return fixed_stages

    include_stages = _rename_input_stages(include_stages)
    exclude_stages = _rename_input_stages(exclude_stages)

    selected_stages = (include_stages or default_enabled) - exclude_stages
    debug(f'Selected_stages: {selected_stages}')
    return selected_stages


def recursively_search_input(input_directory):
    result = {'bcbio': list(), 'dragen': list()}
    directories = [input_directory]
    i = 0
    while directories:
        i += 1
        if i > 1000:
            msg_p1 = f'reached maximum allowed iterations for {input_directory} while'
            msg_p2 = 'searching for input, check that input directories are correct'
            critical(f'ERROR: {msg_p1} {msg_p2}')
        directory = directories.pop()
        if is_bcbio_directory(directory):
            result['bcbio'].append(directory)
        elif is_dragen_directory(directory):
            result['dragen'].append(directory)
        elif directory.is_dir():
            directories.extend(p for p in directory.iterdir() if p.is_dir())
        else:
            assert False
    return result


def is_dragen_directory(path):
    # DRAGEN output directories are identified by the presence of a *-replay.json file
    replay_fps = list(path.glob('*replay.json'))
    return bool(replay_fps)


def is_dragen_tumor_normal_directory(path):
    # DRAGEN tumor/normal directories are differentiated from normal directories by presence of a
    # tumor BAM file with the filename pattern '*_tumor.bam'
    tumor_bam = list(path.glob('*_tumor.bam'))
    return bool(tumor_bam)


def is_bcbio_directory(path):
    # bcbio output directories are identified by the presence of both a 'final' and 'config'
    # directory
    final_dir = path / 'final'
    config_dir = path / 'config'
    return final_dir.is_dir() and config_dir.is_dir()


def get_dragen_output_prefix(dirpath):
    for fp in dirpath.iterdir():
        if not fp.match('*replay.json'):
            continue
        return fp.name.replace('-replay.json', '')
    else:
        critical('could not determine output prefix for DRAGEN directory \'{dirpath}\'')


def create_dragen_paired_directories_from_config(smconfig):
    # Set subject identifier
    tumor_subject_id_inferred = get_subject_id_from_dragen_dir(smconfig['dragen_somatic_dir'])
    normal_subject_id_inferred = get_subject_id_from_dragen_dir(smconfig['dragen_germline_dir'])
    if smconfig.get('dragen_subject_id'):
        subject_id = smconfig.get('dragen_subject_id')
    elif tumor_subject_id_inferred == None:
        critical('could not infer subject id from somatic dir, please specify with --dragen_subject_id')
    elif normal_subject_id_inferred == None:
        critical('could not infer subject id from germline dir, please specify with --dragen_subject_id')
    elif tumor_subject_id_inferred != normal_subject_id_inferred:
        critical(
            f'got different subject ids from somatic ({tumor_subject_id_inferred}) and germline'
            f' ({normal_subject_id_inferred}) dirs, please specify with --dragen_subject_id'
        )
    else:
        subject_id = tumor_subject_id_inferred
    # Set tumor identifier
    tumor_samples_inferred = get_samples_from_dragen_dir_bams(smconfig['dragen_somatic_dir'])
    if smconfig.get('dragen_tumor_id'):
        tumor_id = smconfig.get('dragen_tumor_id')
        if tumor_id != tumor_samples_inferred['tumor']:
            warn(
                f'provided DRAGEN tumor id ({tumor_id}) doesn\'t match id collected'
                f' from discovered BAM file ({tumor_samples_inferred["tumor"]})'
            )
    else:
        tumor_id = tumor_samples_inferred['tumor']
    # Set normal identifier
    normal_samples_inferred = get_samples_from_dragen_dir_bams(smconfig['dragen_germline_dir'])
    if smconfig.get('dragen_normal_id'):
        normal_id = smconfig.get('dragen_normal_id')
        if normal_id != normal_samples_inferred['normal']:
            warn(
                f'provided DRAGEN normal id ({normal_id}) doesn\'t match id collected'
                f' from discovered BAM file ({normal_samples_inferred["normal"]})'
            )
        if 'normal' in tumor_samples_inferred and normal_id != tumor_samples_inferred['normal']:
            warn(
                f'provided DRAGEN normal id ({normal_id}) doesn\'t match id collected'
                f' from discovered BAM file ({tumor_samples_inferred["normal"]})'
            )
    else:
        normal_id = normal_samples_inferred['normal']
    # Create datastructure used for DragenProject init
    return [
        {
            'subject_id': subject_id,
            'tumor_normal_run': {
                'normal': normal_id,
                'tumor': tumor_id,
                'path': smconfig['dragen_somatic_dir'],
                'prefix': get_dragen_output_prefix(smconfig['dragen_somatic_dir'])
            },
            'normal_run': {
                'normal': normal_id,
                'path': smconfig['dragen_germline_dir'],
                'prefix': get_dragen_output_prefix(smconfig['dragen_germline_dir'])
            },
        }
    ]


def pair_dragen_directories(paths):
    # DRAGEN tumor/normal and normal directories are paired on the basis of the normal sample
    # name.
    #
    # Tumor and normal sample names are extracted from the BAM header. Specifically, the BAM
    # sample name is retrieved from the 'SM' (sample) field of the '@RG' (read group) header line.
    #
    # Tumor or normal identity of a sample is inferred from the BAM filename: if a BAM filename
    # contains the '_tumor.bam' suffix then it and the sample name is set as the tumor, otherwise
    # set as the normal sample.
    #
    # The subject identifier is from the DRAGEN output directory name.
    #
    # Assumes a one-to-one pairing for DRAGEN tumor/normal and normal output directories i.e. no
    # multiple tumor/normal runs to a single normal run.

    # Sort paths by normal sample name so that normal and tumor/normal are placed together
    paths_sorted = dict()
    for path in paths:
        dir_type = 'tumor_normal_run' if is_dragen_tumor_normal_directory(path) else 'normal_run'
        samples = get_samples_from_dragen_dir_bams(path)
        # Ensure we have found normal names
        if 'normal' not in samples:
            critical(f'Could not find normal sample name for DRAGEN directory {path}')
        # Sort by normal sample name, add path, subject ID to stored data
        sample_normal = samples['normal']
        if sample_normal not in paths_sorted:
            paths_sorted[sample_normal] = dict()
        assert dir_type not in paths_sorted[sample_normal]
        paths_sorted[sample_normal][dir_type] = samples
        paths_sorted[sample_normal][dir_type]['path'] = path
        paths_sorted[sample_normal][dir_type]['prefix'] = get_dragen_output_prefix(path)
        paths_sorted[sample_normal][dir_type]['subject_id'] = get_subject_id_from_dragen_dir(path)

    # Differentiated paired and unpaired paths
    paths_unpaired = list()
    paths_paired = list()
    for paths in paths_sorted.values():
        if 'normal_run' in paths and 'tumor_normal_run' in paths:
            # Ensure we have collected only one subject id for this set of inputs
            assert len({d['subject_id'] for d in paths.values()}) == 1
            paths['subject_id'] = paths['normal_run']['subject_id']
            paths_paired.append(paths)
        else:
            for dir_type, data in paths.items():
                paths_unpaired.append((dir_type, data['path']))
    # Emit warning for unpaired paths
    if paths_unpaired:
        paths_unpaired_strs = list()
        for dir_type, path in paths_unpaired:
            paths_unpaired_strs.append(f'{dir_type}: {path}')
        paths_unpaired_str = '\n\t'.join(paths_unpaired_strs)
        warn(f'could not pair DRAGEN directories:\n\t{paths_unpaired_str}')
    return paths_paired


def get_samples_from_dragen_dir_bams(dir_fp):
    samples = dict()
    for bam_fp in dir_fp.glob('*bam'):
        # Set BAM type on basis of suffix
        bam_type = 'tumor' if bam_fp.name.endswith('_tumor.bam') else 'normal'
        assert bam_type not in samples
        samples[bam_type] = get_read_group_sample_name(bam_fp)
    return samples


def get_read_group_sample_name(bam_fp):
    bam = pysam.AlignmentFile(bam_fp)
    header = bam.header.to_dict()
    samples = {rg['SM'] for rg in header.get('RG', list())}
    if len(samples) == 0:
        critical(f'could not retrieve sample name from the @RG SM field in {bam_fp}')
    elif len(samples) > 1:
        critical(
            'found more than one sample name in the @RG SM fields for '
            f'{bam_fp}: {", ".join(samples)}'
        )
    return samples.pop()


def get_subject_id_from_dragen_dir(dir_fp):
    # NOTE(SW): this is extremely fragile and will only work with the current naming system. This
    # should be a focus for improvement to increase stability.
    return re.sub('_.+$', '', dir_fp.name)


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


def cnt_vars(vcf_path, passed=False):
    import cyvcf2
    snps = 0
    indels = 0
    others = 0
    for rec in cyvcf2.VCF(vcf_path):
        if passed and rec.FILTER is not None and rec.FILTER != 'PASS':
            continue
        if rec.is_snp:
            snps += 1
        elif rec.is_indel:
            indels += 1
        else:
            others += 1
    return snps, indels, others

def pierian_subset_snvs_cmd(input_vcf, MAX=50_000):
    print("Preparing Pierian SNV subset (if required)")
    def _count_vars(vcf_path, bcftools_filter_expr=None, bed_regions=None, msg=None):
        cmd = f'bcftools view {vcf_path} '
        if bed_regions:
            cmd += f'-R {bed_regions} '
        if bcftools_filter_expr:
            cmd += f' | bcftools filter {bcftools_filter_expr} '
        cmd += ' | bcftools view -H | wc -l'
        ret = int(subprocess.check_output(cmd, shell=True).strip())
        print(f'{msg}{ret}')
        return ret

    expr = {
            'e1': {
                'exp': '-e "(PCGR_TIER == \'NONCODING\') && (PCGR_CONSEQUENCE == \'intergenic_variant\')"',
                'msg': 'dump NONCODING, intergenic_variant: '
                },
            'e2': {
                'exp': '-e "(PCGR_TIER == \'NONCODING\') && (PCGR_CONSEQUENCE == \'intergenic_variant\' || PCGR_CONSEQUENCE == \'intron_variant\')"',
                'msg':  'dump NONCODING, intergenic_variant or intron_variant: '
                },
            'e3': {
                'exp': '-e "(PCGR_TIER == \'NONCODING\') && (PCGR_CONSEQUENCE == \'intergenic_variant\' || PCGR_CONSEQUENCE == \'intron_variant\' || PCGR_CONSEQUENCE == \'intron_variant|_non_coding_transcript_variant\')"',
                'msg':  'dump NONCODING, intergenic_variant or intron_variant or intron_variant|_non_coding_transcript_variant: '
                }
            }

    cmd = ''
    msg = ''
    bed = get_all_genes_bed()
    if _count_vars(input_vcf, msg='Total: ') < MAX:
        msg = "No Pierian SNV subsetting required!"
        cmd = f'bcftools view {input_vcf}'
        return({'msg': msg, 'cmd': cmd})
    else:
        # subset to gene regions
        cmd = f'bcftools view -R {bed} {input_vcf} | '

    if _count_vars(input_vcf, bed_regions=bed, msg='Total in gene regions: ') < MAX:
        msg = "just subsetting to gene regions"
        cmd += f'bcftools view'
    elif _count_vars(input_vcf, bed_regions=bed, bcftools_filter_expr=expr["e1"]["exp"], msg=expr["e1"]["msg"]) < MAX:
        msg = "do not include noncoding intergenic"
        cmd += f'bcftools filter {expr["e1"]["exp"]}'
    elif _count_vars(input_vcf, bed_regions=bed, bcftools_filter_expr=expr["e2"]["exp"], msg=expr["e2"]["msg"]) < MAX:
        msg = "do not include noncoding intergenic/intronic"
        cmd += f'bcftools filter {expr["e2"]["exp"]}'
    elif _count_vars(input_vcf, bed_regions=bed, bcftools_filter_expr=expr["e3"]["exp"], msg=expr["e3"]["msg"]) < MAX:
        msg = "do not include noncoding intergenic/intronic/_non_coding_transcript_variant"
        cmd += f'bcftools filter {expr["e3"]["exp"]}'
    else:
        msg = "just filter out all noncoding"
        cmd += f'bcftools filter -e "PCGR_TIER == \'NONCODING\'"'
    return({'msg': msg, 'cmd': cmd})
