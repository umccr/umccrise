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
    Loc = collections.namedtuple('Loc', 'name host_pattern hsapiens extras')
    hostname = socket.gethostname()
    for loc in [
        Loc(name='spartan',
            host_pattern=r'spartan.*\.hpc\.unimelb\.edu\.au',
            hsapiens='/home/vlad/bcbio/genomes/Hsapiens',
            extras='/data/projects/punim0010/local/share/extras/',
        ),
        Loc(name='raijin',
            host_pattern=r'^raijin',
            hsapiens='/home/563/vs2870/g/bcbio/genomes/Hsapiens',
            extras='/g/data3/gx8/extras/',
        ),
        Loc(name='vlad',
            host_pattern=r'^5180L-135800-M.local$',
            hsapiens='/Users/vsaveliev/genomes/Hsapiens',
            extras='/Users/vsaveliev/Analysis/umccrize/',
        ),]:
        if re.match(loc.host_pattern, hostname):
            return loc
    raise Exception('Could not find loc for hostname ' + hostname)


def get_sample_name(vcf_path):
    with open_gzipsafe(vcf_path) as f:
        for line in f:
            m = re.match(r'^##SAMPLE=<ID=(?P<name>\S+),Genomes=Tumor>$', line)
            if m:
                return m.group('name')
            m = re.match(r'^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t(?P<name>\S+)', line)
            if m:
                return m.group('name')
    raise ValueError


BCINSTALL = "/data/projects/punim0010/local/share/bcbio/"
BCRESULT = "/data/cephfs/punim0010/data/Results/Avner/WPT-013/final/"
BCFINAL = "2017-10-19_WPT-013"
BCPOST = "/data/cephfs/punim0010/projects/Hofmann_Explore/testrun"
EXTRAS = "/data/projects/punim0010/local/share/extras/"

# This information should come from the *.csv
BCTUMOR = "WPT-013-organoid"
BCNORMAL = "WPT-013-normal"
BCBATCH = "batch1"
