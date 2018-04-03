import os
import sys
from os.path import isfile, join, dirname, abspath
from ngs_utils.file_utils import verify_file
import click
import subprocess

import locale
try:
    if 'UTF-8' not in locale.getlocale(locale.LC_ALL):
        locale.setlocale(locale.LC_ALL, 'en_AU.UTF-8')
except TypeError:
    pass


def package_path():
    return dirname(abspath(__file__))


@click.command()
@click.argument('bcbio_project', type=click.Path(exists=True))
@click.argument('rule', nargs=-1)
@click.option('-o', 'output_dir', type=click.Path())
@click.option('-j', 'jobs', default=1)
@click.option('-s', '--sample', 'sample')
@click.option('-b', '--batch', 'batch')
@click.option('-u', '--uid', '--uuid', 'unique_id', default='XXXXXX')
@click.option('--unlock', is_flag=True)
def main(bcbio_project, rule=list(), output_dir=None, jobs=1, sample=None, batch=None, unique_id=None, unlock=False):
    rule = list(rule)

    bcbio_project = os.path.abspath(bcbio_project)

    conf = f'run_dir={bcbio_project}'

    if sample:
        conf += f' sample={sample}'
    if batch:
        conf += f' batch={batch}'

    # if 'pcgr_download' in rule and not unique_id:
    #     sys.stderr.write(f'Error: when you run pcgr_download, provide the unique id with --uid option so umccrise can find the tarballs:\n')
    #     sys.stderr.write('\n')
    #     args = ' '.join(sys.argv)
    #     sys.stderr.write(f'    {args} --uid XXXXXX\n')
    #     sys.stderr.write('\n')
    #     sys.exit(1)
    # if unique_id:
    conf += f' unique_id="{unique_id}"'
    conf += f' pcgr_download=yes'

    output_dir = output_dir or 'umccrised'
    output_dir = abspath(output_dir)

    cmd = (f'snakemake ' +
        f'{" ".join(rule)} ' +
        f'--snakefile {join(package_path(), "Snakefile")} ' +
        f'--printshellcmds ' +
        f'--directory {output_dir} ' +
        f'-j {jobs} ' +
        f'--config {conf} ')

    if unlock:
        print('* Unlocking previous run... *')
        print(cmd + ' --unlock')
        subprocess.call(cmd + ' --unlock', shell=True)
        print('* Now rerunning *')

    print(cmd)
    exit_code = subprocess.call(cmd, shell=True)
    if exit_code != 0:
        sys.stderr.write('--------\n')
        sys.stderr.write(f'Error running Umccrise: snakemake returned a non-zero status.\n')
        sys.exit(exit_code)

    # Cleanup
    work_dir = join(output_dir, 'work')
    # if isdir(work_dir):
    #     shutils.rmtree(work_dir)


def get_sig_rmd_file():
    """ Returns path to sig.Rmd file - R-markdown source for mutational signature analysys.
        The file must be located at the same directory as the Snakefile and the patient_analysis module.
    """
    return verify_file(join(package_path(), 'sig.Rmd'))

def get_signatures_probabilities():
    return verify_file(join(package_path(), 'rmd_files', 'signatures_probabilities.txt'))

def get_suppressors():
    return verify_file(join(package_path(), 'rmd_files', 'suppressors.txt'))

def get_cancer_genes_ensg():
    return verify_file(join(package_path(), 'cancer_genes_ENSG.txt'))
