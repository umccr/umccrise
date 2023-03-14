#!/usr/bin/env python
from setuptools import setup, find_packages
from versionpy import get_version, find_package_files

pkg = 'umccrise'

version = get_version(pkg)

setup(
    name='umccrise',
    version=version,
    author='Vlad Saveliev',
    author_email='vladislav.sav@gmail.com',
    description='UMCCRisation of DRAGEN analysis results',
    keywords='bioinformatics',
    url='https://github.com/umccr/umccrise',
    license='GPLv3',
    package_data={
        pkg: find_package_files('', pkg),
    },
    packages=find_packages(),
    scripts=[
        'scripts/umccrise',
        'scripts/pcgr_wrap',
        'scripts/cpsr_wrap',
        'scripts/conpair',
        'scripts/umccrise_refdata_pull',
    ],
    include_package_data=True,
    zip_safe=False,
    # For MultiQC_umccr
    entry_points = {
        'multiqc.templates.v1': [
            'umccr = umccrise.multiqc.templates.umccr',
        ],
        'multiqc.hooks.v1': [
            'config_loaded            = umccrise.multiqc.multiqc_umccr:config_loaded',
            'execution_start          = umccrise.multiqc.multiqc_umccr:execution_start',
            'after_modules            = umccrise.multiqc.multiqc_umccr:before_set_general_stats_html',
            'before_report_generation = umccrise.multiqc.multiqc_umccr:after_set_general_stats_html',
        ]
    },
)
