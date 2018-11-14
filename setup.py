#!/usr/bin/env python
from setuptools import setup, find_packages
from releazit import get_version, find_package_files

pkg = 'umccrise'

version = get_version(pkg)

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
        pkg: find_package_files('', pkg),
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
