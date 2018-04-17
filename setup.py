#!/usr/bin/env python
import sys
import os
from os.path import join, isfile, abspath, dirname
from setuptools import setup

version = '0.4.0'

setup(
    name='umccrise',
    version=version,
    author='Vlad Saveliev',
    description='UMCCRisation of bcbio-nextgen analysis results',
    keywords='bioinformatics',
    license='GPLv3',
    packages=[
        'umccrise',
    ],
    scripts=[
        'scripts/vcfToBedpe',
    ],
    entry_points = {
        'console_scripts': [
            'umccrise=umccrise:main',
            'umccrize=umccrise:main',
        ],
    },
    include_package_data=True,
)
