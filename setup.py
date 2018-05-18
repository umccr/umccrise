#!/usr/bin/env python
from setuptools import setup

version = '0.6.0'

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
    ],
    entry_points = {
        'console_scripts': [
            'umccrise=umccrise:main',
            'umccrize=umccrise:main',
            'pcgr=umccrise.pcgr.pcgr_runner:main',
        ],
    },
    include_package_data=True,
)
