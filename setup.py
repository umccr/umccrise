#!/usr/bin/env python
from setuptools import setup, find_packages
from ngs_utils import setup_utils

name = 'umccrise'

version = setup_utils.get_cur_version(name)

setup(
    name='umccrise',
    version=version,
    author='Vlad Saveliev',
    author_email='vladislav.sav@gmail.com',
    description='UMCCRisation of bcbio-nextgen analysis results',
    keywords='bioinformatics',
    url='https://github.com/umccr/umccrise',
    license='GPLv3',
    package_data={
        name: setup_utils.find_package_files('', name)
    },
    packages=find_packages(),
    scripts=[
        'vendor/vcfToBedpe',
        'scripts/umccrise',
        'scripts/pcgr',
    ],
    include_package_data=True,
    zip_safe=False,
)
