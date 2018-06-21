#!/usr/bin/env python
from setuptools import setup

version = '0.7'

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
