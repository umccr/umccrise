#!/usr/bin/env python
import sys
import os
from os.path import join, isfile, abspath, dirname
from setuptools import setup

version = '0.1'

setup(
    name='somatic_filtering',
    version=version,
    author='Vlad Saveliev',
    description='VCF file normalization, validation and filtering',
    keywords='bioinformatics',
    license='GPLv3',
    packages=[
        'somatic_filtering'
    ],
    scripts=[
        'scripts/normalize.py',
        'somatic_filtering/vardict/proc_vardict_vcf'
    ],
    include_package_data=True,
    install_requires=[],
)
