#!/usr/bin/env python

import sys
from pprint import pprint

import click
from os.path import isfile, join
from cyvcf2 import VCF, Writer
from ngs_utils.file_utils import verify_file
from ngs_utils.vcf_utils import get_sample_ids


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-o', 'output_file', type=click.Path())
def main(input_file, output_file=None):
    vcf = VCF(input_file, gts012=True)

    if output_file:
        w = Writer(output_file, vcf)
    else:
        w = None
        sys.stdout.write('\n'.join([h for h in vcf.raw_header.split('\n') if not h.startswith('##SAMPLE')]))

    tumor_index, control_index = get_sample_ids(input_file)

    header_by_tag = _find_tags(vcf)

    # Add headers
    for t, h in header_by_tag.items():
        t = t.upper()
        vcf.add_info_to_header({
            'ID': t,
            'Description': f'{t} in tumor sample',
            'Type': 'Float' if t == 'AF' else 'Integer',
            'Number': '1'
        })
        if control_index:
            t = 'NORMAL_' + t
            vcf.add_info_to_header({
                'ID': t,
                'Description': f'{t} in control sample',
                'Type': 'Float' if t == 'AF' else 'Integer',
                'Number': '1'
            })

    # Go through each record and add new INFO fields
    for rec in vcf:
        val_per_tag_per_sample = _collect_vals_per_sample(rec, control_index, tumor_index, header_by_tag)

        for tag, val_by_sample in val_per_tag_per_sample.items():
            tag = tag.upper()
            rec.INFO[tag] = val_by_sample['tumor']
            if val_by_sample.get('normal'):
                rec.INFO['NORMAL_' + tag] = val_by_sample['normal']

        if w:
            w.write_record(rec)
        else:
            sys.stdout.write(str(rec))

    if w:
        w.close()
    vcf.close()


def _find_tags(vcf):
    tags_options = {
        'ad': {'AD'},  # contains raw counts for each allele
        'af': {'AF', 'FREQ', 'FA'},  # contains actual frequencies which can be used directly
        'dp': {'DP'},
        'mq': {'MQ', 'MMQ'},
    }

    formats_by_id = {}
    infos_by_id = {}
    for h in vcf.header_iter():
        i = h.info()
        ht = i['HeaderType']
        if ht in ['FORMAT', 'INFO']:
            if i['Type'] == 'Integer':
                i['python_type'] = int
            if i['Type'] == 'Float':
                i['python_type'] = float

            if i['HeaderType'] == 'FORMAT':
                formats_by_id[i['ID']] = i
            elif i['HeaderType'] == 'INFO':
                infos_by_id[i['ID']] = i

    header_by_tag = {t: None for t in tags_options}

    for headers_by_id in [formats_by_id, infos_by_id]:
        for tag, header in header_by_tag.items():
            if header is None:
                header_id_options = tags_options[tag].intersection(headers_by_id)
                if header_id_options:
                    header_by_tag[tag] = headers_by_id[header_id_options.pop()]

    pprint(header_by_tag)
    return header_by_tag


def _collect_vals_per_sample(rec, control_index, tumor_index, header_by_tag):
    d = {t: {'tumor': -1, 'normal': -1} for t in header_by_tag}  # val_by_sample_by_tags

    if header_by_tag['af']:
        af_header = header_by_tag['af']
        ad_tag_is_frequency = True
    else:
        af_header = header_by_tag['ad']
        ad_tag_is_frequency = False

    c_alt_support = -1
    t_alt_support = -1
    if af_header:
        if af_header['HeaderType'] == 'FORMAT':
            af_data = rec.format(af_header['ID'])
            t_alt_support = str(af_data[tumor_index][0])

            sample_dim = af_data.shape[0]
            if sample_dim >= 2:
                c_alt_support = str(af_data[control_index][0])

        else:  # Strelka doesn't have AF or AD in FORMAT, using INFO["AF"] (unfortunately INFO["AF"] is always 0.25)
            t_alt_support = str(rec.INFO[af_header['ID']])

    if t_alt_support == '-2147483648': t_alt_support = -1
    if c_alt_support == '-2147483648': c_alt_support = -1

    dp_header = header_by_tag['dp']
    if dp_header:
        if dp_header['HeaderType'] == 'FORMAT':
            dp_data = rec.format(dp_header['ID'])
            d['dp']['tumor'] = str(dp_data[tumor_index][0])

            sample_dim = dp_data.shape[0]
            if sample_dim >= 2:
                d['dp']['normal'] = str(dp_data[control_index][0])

            if d['dp']['tumor']  == '-2147483648': d['dp']['tumor']  = -1
            if d['dp']['normal'] == '-2147483648': d['dp']['normal'] = -1
        else:
            d['dp']['tumor'] = str(rec.INFO[dp_header['ID']])

    if ad_tag_is_frequency:
        d['af']['tumor'] = t_alt_support
    elif t_alt_support != -1 and d['dp']['tumor'] > 0:
        d['af']['tumor'] = '{0:.4f}'.format(int(t_alt_support) / float(d['dp']['tumor']))

    if ad_tag_is_frequency:
        d['af']['normal'] = c_alt_support
    elif c_alt_support != -1 and d['dp']['normal'] != -1:
        if float(d['dp']['normal']) > 0:
            d['af']['normal'] = "{0:.4f}".format(int(c_alt_support) / float(d['dp']['normal']))

    _parse_tag(rec, header_by_tag, 'mq', d, tumor_index, control_index)

    return d


def _parse_tag(rec, header_by_tag, tag, d, tumor_index, control_index):
    header = header_by_tag[tag]
    if header:
        if header['HeaderType'] == 'FORMAT':
            data = rec.format(header['ID'], header['python_type'])
            d[tag]['tumor'] = str(data[tumor_index][0])

            sample_dim = data.shape[0]
            if sample_dim >= 2:
                d[tag]['normal'] = str(data[control_index][0])

            if d[tag]['tumor']  == '-2147483648': d[tag]['tumor']  = -1
            if d[tag]['normal'] == '-2147483648': d[tag]['normal'] = -1
        else:
            d[tag]['tumor'] = str(rec.INFO[header['ID']])


if __name__ == '__main__':
    main()
