#!/usr/bin/env python
from setuptools import setup

with open('VERSION.txt') as f:
    version = f.read().strip().split('\n')[0]

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
        'vendor/vcfToBedpe',
        'scripts/umccrise',
        'scripts/pcgr',
    ],
    include_package_data=True,
)
