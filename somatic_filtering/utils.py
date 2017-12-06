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


def get_ref_path():
    hostname = socket.gethostname()
    if hostname == '5180L-135800-M.local':
        loc = 'local'
        ref_loc = '/Users/vsaveliev/genomes/Hsapiens'
    elif re.match(r'spartan.*\.hpc\.unimelb\.edu\.au', hostname):
        loc = 'spartan'
        ref_loc = '/home/vlad/bcbio/genomes/Hsapiens'
    elif hostname.startswith('raijin'):
        loc = 'raijin'
        ref_loc = '/home/563/vs2870/g/bcbio/genomes/Hsapiens'
    else:
        ref_loc = ''
    return ref_loc


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


