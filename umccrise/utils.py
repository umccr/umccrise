import collections
import contextlib
import csv
import hashlib
import os
import shutil
import socket
import subprocess
import time
import re
import gzip
import yaml
from ngs_utils.file_utils import open_gzipsafe


##############################
### HPC dependencies paths ###


def find_loc():
    """ Depending on the machine name, return a dict conatining system-dependant paths
        to human reference genomes and extras
    """
    Loc = collections.namedtuple('Loc',
         'name '
         'host_pattern '
         'hsapiens '
         'extras '
         'panel_of_normals_dir '
         'truth_sets ')

    hostname = socket.gethostname()
    loc_by_name = {
        'spartan': Loc(
            name='spartan',
            host_pattern=r'spartan.*\.hpc\.unimelb\.edu\.au',
            hsapiens='/data/projects/punim0010/local/share/bcbio/genomes/Hsapiens',
            extras='/data/cephfs/punim0010/extras',
            panel_of_normals_dir='/data/cephfs/punim0010/extras/panel_of_normals',
            truth_sets={
                'mb': {
                    'GRCh37': {
                        'vcf': '/data/cephfs/punim0010/data/External/Reference/ICGC_MB/MB-benchmark.vcf.gz',
                    }
                },
                'dream': {
                    'GRCh37': {
                        'vcf': 'GRCh37/validation/dream-syn3/truth_small_variants.vcf.gz',
                        # 'bed': 'GRCh37/validation/dream-syn3/truth_regions.bed',
                    }
                },
                'giab': {
                    'GRCh37': {
                        'vcf': 'GRCh37/validation/giab-NA12878/truth_small_variants.vcf.gz',
                        'bed': 'GRCh37/validation/giab-NA12878/truth_regions.bed',
                    }
                },
                'colo': {
                    'GRCh37': {
                        'vcf': '/data/cephfs/punim0010/data/External/Reference/COLO829_Craig/truth_set/EGAZ00001226241_ListforNatureReports.IndelsandSNVs.final.Suppl1.snpEff.validated.SORTED.vcf'
                    }
                }
            },
        ),
        'raijin': Loc(
            name='raijin',
            host_pattern=r'^raijin|(r\d\d\d\d$)',
            hsapiens='/g/data/gx8/local/development/bcbio/genomes/Hsapiens',
            extras='/g/data3/gx8/extras',
            panel_of_normals_dir='/g/data3/gx8/extras/panel_of_normals',
            truth_sets={
                'giab': {
                    'GRCh37': {
                        'vcf': 'GRCh37/validation/giab-NA12878/truth_small_variants.vcf.gz',
                        'bed': 'GRCh37/validation/giab-NA12878/truth_regions.bed',
                    }
                }
            }
        ),
        'vlad': Loc(
            name='vlad',
            host_pattern=r'^5180L-135800-M.local$',
            hsapiens='/Users/vsaveliev/genomes/Hsapiens',
            extras='/Users/vsaveliev/Analysis/umccrise',
            panel_of_normals_dir='/Users/vsaveliev/Analysis/panel_of_normals/GRCh37/normals',
            truth_sets={
                'giab': {
                    'GRCh37': {
                        'bed': 'GRCh37/validation/giab-NA12878/truth_regions.bed',
                    }
                }
            }
        ),
        'travis': Loc(
            name='travis',
            host_pattern=r'^travis-',
            hsapiens='../../data/genomes/Hsapiens',
            # For tests, using the following in Hsapiens:
            # For germline subsampling: 'GRCh37/coverage/prioritize/cancer/az300.bed.gz'
            # For goleft depth:         'GRCh37/seq/GRCh37.fa
            # For bedtools slop :       'GRCh37/seq/GRCh37.fa.fai'
            extras='',
            panel_of_normals_dir='../../data/panel_of_normals',
            truth_sets={
                'giab': {
                    'GRCh37': {
                        'bed': 'GRCh37/validation/giab-NA12878/truth_regions.bed',
                    }
                }
            }
        ),
    }
    if 'TRAVIS' in os.environ.keys():
        return loc_by_name['travis']
    else:
        for loc in loc_by_name.values():
            if re.match(loc.host_pattern, hostname):
                return loc

    return None


def get_loc():
    loc = find_loc()
    if loc:
        return loc
    else:
        raise Exception('Could not find loc for hostname ' + socket.gethostname())


"""
All known truth sets
"""
truth_vcfs = {
    'spartan': {
        'mb': {
            'GRCh37': {
                'vcf': '/data/cephfs/punim0010/data/External/Reference/ICGC_MB/MB-benchmark.vcf.gz'
            }
        },
        'dream': {
            'GRCh37': {
                'vcf': 'GRCh37/validation/dream-syn3/truth_small_variants.vcf.gz',
                # 'bed': 'GRCh37/validation/dream-syn3/truth_regions.bed'
            }
        }
    }
}


#################
### VCF utils ###

def get_sample_names(vcf_path):
    return get_sample_ids(vcf_path)

def get_tumor_sample_name(vcf_path):
    return get_sample_ids(vcf_path, return_names=True)[0]

def get_normal_sample_name(vcf_path):
    return get_sample_ids(vcf_path, return_names=True)[1]

def get_tumor_sample_id(vcf_path):
    return get_sample_ids(vcf_path)[0]

def get_normal_sample_id(vcf_path):
    return get_sample_ids(vcf_path)[1]

def get_sample_ids(vcf_path, return_names=False):
    """ Finds tumor and control sample names/ids from a bcbio-derived VCF,
        and returns a tuple (tumor, control)

        1. If ##SAMPLE header field is found, use that for the name.
        2. Otherwise, return the first (for tumor) and second (for normal)
           sample in #CHROM header.
    """
    tumor_name, control_name = None, None

    with open_gzipsafe(vcf_path) as f:
        for line in f:
            m = re.match(r'^##SAMPLE=<ID=(?P<name>\S+),Genomes=Tumor>$', line)
            if m:
                tumor_name = m.group('name')
            m = re.match(r'^##SAMPLE=<ID=(?P<name>\S+),Genomes=Germline>$', line)
            if m:
                control_name = m.group('name')

            if tumor_name and return_names:
                return tumor_name, control_name

            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                tumor_id, control_id = None, None
                if tumor_name:
                    tumor_id = samples.index(tumor_name)
                if control_name:
                    control_id = samples.index(control_name)

                if tumor_id is None and control_id is not None:
                    tumor_id = 1 if control_id == 0 else 0

                elif control_id is None and tumor_id is not None and len(samples) > 1:
                    control_id = 1 if tumor_id == 0 else 0

                elif control_id is None and tumor_id is None:
                    tumor_id = 0
                    if len(samples) > 1:
                        control_id = 1

                if tumor_name is None:
                    tumor_name = samples[tumor_id]
                    if control_name is None and len(samples) > 1:
                        control_name = samples[control_id]

                if return_names:
                    return tumor_name, control_name
                else:
                    return tumor_id, control_id
    raise ValueError
