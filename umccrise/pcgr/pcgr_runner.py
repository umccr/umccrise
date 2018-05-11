import sys
from os import rename
from os.path import isfile, join, dirname, abspath, basename
from ngs_utils.file_utils import verify_file, safe_mkdir, splitext_plus
from ngs_utils import logger
import click
import subprocess
from python_utils.hpc import find_loc
from umccrise import package_path


@click.command()
@click.argument('vcf_path', type=click.Path(exists=True))
@click.argument('cnv_path', type=click.Path(exists=True), required=False)
@click.option('-g', 'genome')
@click.option('-o', 'output_dir', type=click.Path())
@click.option('-s', 'sample')
@click.option('--germline', is_flag=True)
@click.option('--no-docker', is_flag=True)
def main(vcf_path, cnv_path=None, output_dir=None, genome='GRCh37', sample=None, germline=False, no_docker=False):

    loc = find_loc()
    pcgr_dir = loc.pcgr_dir
    if not pcgr_dir:
        logger.critical(f'PCGR is not installed on the system "{loc.name}".'
                        f' Please, pull it with `git clone https://github.com/vladsaveliev/pcgr`'
                        f' and install with `bash pcgr/install_no_docker/install.sh`')
    output_dir = output_dir or 'pcgrred'
    output_dir = abspath(output_dir)
    safe_mkdir(output_dir)

    somatic_toml = join(package_path(), 'pcgr', 'pcgr_configuration_somatic.toml')
    germline_toml = join(package_path(), 'pcgr', 'pcgr_configuration_normal.toml')

    logger.init(log_fpath_=join(output_dir, 'pcgr.log'), save_previous=True)

    sample = sample or splitext_plus(basename(vcf_path))[0]
    pcgr_genome = "grch38" if genome in ["hg38", "GRCh38"] else "grch37"
    expected_report_path = join(output_dir, f'{sample}.pcgr_acmg.{pcgr_genome}.html')
    renamed_report_path = join(output_dir, f'{sample}.pcgr_acmg.html')

    cmd = (f'{join(pcgr_dir, "pcgr.py")}'
           f' --input_vcf {abspath(vcf_path)}'
           f' {("--input_cna " + abspath(cnv_path)) if cnv_path else ""}'
           f' {pcgr_dir}'
           f' {output_dir}'
           f' {pcgr_genome}'
           f' {somatic_toml if not germline else germline_toml}'
           f' {sample}'
           f' {" --no-docker " if no_docker else " --docker-uid root"}'
           f' --force_overwrite'
    )
    if isfile(expected_report_path):
        rename(expected_report_path, renamed_report_path)

    print(cmd)
    exit_code = subprocess.call(cmd, shell=True)
    if exit_code != 0:
        sys.stderr.write('--------\n')
        sys.stderr.write(f'Error running PCGR.\n')
        sys.exit(exit_code)
