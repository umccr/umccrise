#!/usr/bin/env python

import sys
import click
from os.path import isfile, join
from cyvcf2 import VCF, Writer
from ngs_utils.file_utils import verify_file
from ngs_utils.vcf_utils import get_tumor_sample_name, get_normal_sample_name


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-o', 'output_file', type=click.Path())
def main(input_file, output_file=None):
    tumor_sample_name = get_tumor_sample_name(input_file)
    control_sample_name = get_normal_sample_name(input_file)
    pull_genotype_data(input_file, control_sample_name, tumor_sample_name, output_file)


TUMOR_AF = 'TVAF'
NORMAL_AF = 'CVAF'
TUMOR_DP = 'TDP'
NORMAL_DP = 'CDP'


def pull_genotype_data(query_vcf, control_sample_name, tumor_sample_name, output_vcf=None):
    vcf = VCF(query_vcf, gts012=True)

    vcf.add_info_to_header({
        'ID': 'CVAF',
        'Description':
        'Allelic fraction of alternate allele in control sample',
        'Type': 'Float',
        'Number': '1'
    })
    vcf.add_info_to_header({
        'ID': 'TVAF',
        'Description': 'Allelic fraction of alternate allele in tumor sample',
        'Type': 'Float',
        'Number': '1'
    })
    vcf.add_info_to_header({
        'ID': 'CDP',
        'Description': 'Depth at variant site in control sample',
        'Type': 'Integer',
        'Number': '1'
    })
    vcf.add_info_to_header({
        'ID': 'TDP',
        'Description': 'Depth at variant site in tumor sample',
        'Type': 'Integer',
        'Number': '1'
    })

    if output_vcf:
        w = Writer(output_vcf, vcf)
    else:
        w = None
        sys.stdout.write(vcf.raw_header)

    try:
        control_index = vcf.samples.index(control_sample_name)
        tumor_index = vcf.samples.index(tumor_sample_name)
    except ValueError as err:
        print("Sample not found in vcf")
        raise

    # Go through each record
    for rec in vcf:
        allelic_support = collect_allelic_support(rec, control_index, tumor_index)

        for tag, val in allelic_support.items():
            if val != -1:
                rec.INFO[tag] = val
            else:
                rec.INFO[tag] = '.'
        if w:
            w.write_record(rec)
        else:
            sys.stdout.write(str(rec))

    if w:
        w.close()
    vcf.close()


def collect_allelic_support(rec, control_index, tumor_index):
    allelic_support = {tag: -1 for tag in [TUMOR_AF, NORMAL_AF, TUMOR_DP, NORMAL_DP]}

    ad_keys = {"AD"}  # contains raw counts
    af_keys = {"AF", "FREQ", "FA"}  # contains actual frequencies which can be used directly
    dp_key = 'DP'

    control_alt_support = -1
    tumor_alt_support = -1

    ad_format_tag = None
    ad_format_tags = af_keys.intersection(rec.FORMAT)
    if ad_format_tags:
        ad_format_tag = ad_format_tags.pop()
        ad_tag_is_frequency = True
    else:
        ad_tag_is_frequency = False
        ad_format_tags = ad_keys.intersection(rec.FORMAT)
        if ad_format_tags:
            ad_format_tag = ad_format_tags.pop()

    if ad_format_tag:
        sample_dat_ao = rec.format(ad_format_tag)
        if sample_dat_ao is not None:
            sample_dim = sample_dat_ao.shape[0]
            if sample_dim >= 2:
                control_alt_support = str(sample_dat_ao[control_index][0])
                if control_alt_support == '-2147483648':
                    control_alt_support = -1
                tumor_alt_support = str(sample_dat_ao[tumor_index][0])
                if tumor_alt_support == '-2147483648':
                    tumor_alt_support = -1
    else:  # Strelka doesn't have AF or AD in FORMAT; reading INFO['AF'] for tumor only instead
        af_keys = af_keys.intersection(rec.INFO)
        if af_keys:
            af_key = af_keys.pop()
            tumor_alt_support = str(rec.INFO[af_key])

    sample_dat_dp = rec.format(dp_key)
    if sample_dat_dp is not None:
        sample_dim = sample_dat_dp.shape[0]
        if sample_dim >= 2:
            allelic_support[NORMAL_DP] = str(sample_dat_dp[control_index][0])
            if allelic_support[NORMAL_DP] == '-2147483648':
                allelic_support[NORMAL_DP] = -1
            allelic_support[TUMOR_DP] = str(sample_dat_dp[tumor_index][0])
            if allelic_support[TUMOR_DP] == '-2147483648':
                allelic_support[TUMOR_DP] = -1

    if ad_tag_is_frequency:
        allelic_support[TUMOR_AF] = tumor_alt_support
    elif tumor_alt_support != -1 and allelic_support[TUMOR_DP] != -1:
        allelic_support[TUMOR_AF] = "{0:.4f}".format(int(tumor_alt_support) / float(allelic_support[TUMOR_DP]))

    if ad_tag_is_frequency:
        allelic_support[NORMAL_AF] = control_alt_support
    elif control_alt_support != -1 and allelic_support[NORMAL_DP] != -1:
        if float(allelic_support[NORMAL_DP]) == 0:  # prevent divide by zero errors
            allelic_support[NORMAL_AF] = 0.0
        else:
            allelic_support[NORMAL_AF] = "{0:.4f}".format(int(control_alt_support) / float(allelic_support[NORMAL_DP]))

    return allelic_support


if __name__ == '__main__':
    main()
