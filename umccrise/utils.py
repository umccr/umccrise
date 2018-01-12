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


def get_loc():
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
    for loc in [
        Loc(name='spartan',
            host_pattern=r'spartan.*\.hpc\.unimelb\.edu\.au',
            hsapiens='/home/vlad/bcbio/genomes/Hsapiens',
            extras='/data/cephfs/punim0010/extras',
            panel_of_normals_dir='/data/cephfs/punim0010/extras/panel_or_normals',
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
                }
            },
        ),
        Loc(name='raijin',
            host_pattern=r'^raijin|(r\d\d\d\d$)',
            hsapiens='/g/data/gx8/local/development/bcbio/genomes/Hsapiens',
            extras='/g/data3/gx8/extras',
            panel_of_normals_dir='/g/data3/gx8/extras/panel_or_normals',
            truth_sets={
                'giab': {
                    'GRCh37': {
                        'vcf': 'GRCh37/validation/giab-NA12878/truth_small_variants.vcf.gz',
                        'bed': 'GRCh37/validation/giab-NA12878/truth_regions.bed',
                    }
                }
            }
        ),
        Loc(name='vlad',
            host_pattern=r'^5180L-135800-M.local$',
            hsapiens='/Users/vsaveliev/genomes/Hsapiens',
            extras='/Users/vsaveliev/Analysis/umccrise',
            panel_of_normals_dir='/Users/vsaveliev/Analysis/panel_of_normals/GRCh37/normals',
            truth_sets={},
        ),]:
        if re.match(loc.host_pattern, hostname):
            return loc
    raise Exception('Could not find loc for hostname ' + hostname)




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

