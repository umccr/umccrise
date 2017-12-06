#!/usr/bin/env python
"""
Normalizes VCF: 
- splits multiallelic ALT
- splits biallelic MNP
- left-aligns indels
- fixes FORMAT and INFO fields
"""
import os
import click
from somatic_filtering.normalization import make_normalize_cmd


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-o', 'output_file', type=click.Path())
@click.option('-f', 'reference_fasta', type=click.Path())
def main(input_file, output_file, reference_fasta=False):
    cmd = make_normalize_cmd(input_file, output_file, reference_fasta)
    print(cmd)
    os.system(cmd)


if __name__ == '__main__':
    main()
