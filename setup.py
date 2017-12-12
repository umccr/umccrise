#!/usr/bin/env python
import sys
import os
from os.path import join, isfile, abspath, dirname
from setuptools import setup

version = '0.1'

setup(
    name='umccrise',
    version=version,
    author='Vlad Saveliev',
    description='VCF file normalization, validation and filtering',
    keywords='bioinformatics',
    license='GPLv3',
    packages=[
        'umccrise',
    ],
    scripts=[
        'scripts/normalize_vcf',
        'scripts/umccrise',
    ],
    include_package_data=True,
)
