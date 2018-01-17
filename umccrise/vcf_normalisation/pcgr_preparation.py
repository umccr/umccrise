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

    # Add headers
    for t in ['AF', 'DP', 'MQ']:
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
        af, dp, mq = _collect_vals_per_sample(rec, control_index, tumor_index)

        for t, v in zip(['AF', 'DP', 'MQ'], [af, dp, mq]):
            rec.INFO[t] = str(v[tumor_index])
            if control_index and v[control_index] is not None:
                rec.INFO['NORMAL_' + t] = str(v[control_index])

        if w:
            w.write_record(rec)
        else:
            sys.stdout.write(str(rec))

    if w:
        w.close()
    vcf.close()


''' 
VarDict
FORMAT/DP,    FORMAT/AF,                           FORMAT/MQ (somatic), INFO/MQ (germline)
             
Mutect2             
FORMAT/DP,    FORMAT/AF,                           FORMAT/MMQ
             
Freebayes                                       
FORMAT/DP     FORMAT/AD = ref_count,alt_count      INFO/MQM+INFO/MQMR
             
GATK-Haplotype                                  
FORMAT/DP     FORMAT/AD = ref_count,alt_count      FORMAT/MMQ
    
Strelka - germline                 
SNV:
FORMAT/DP     FORMAT/AD = ref_count,alt_count      INFO/MQ
INDEL:
FORMAT/DPI    FORMAT/AD = ref_count,alt_count      INFO/MQ

Strelka - somatic                 
SNV:
FORMAT/DP     FORMAT/{ALT}U[0] = alt_count         INFO/MQ
INDEL:
FORMAT/DP     FORMAT/TIR = alt_count               INFO/MQ
'''


def _collect_vals_per_sample(rec, control_index, tumor_index):
    try:
        dp = rec.format('DP')[:,0]
    except:  # Strelka2 indels:
        dp = rec.format('DPI')[:,0]

    # If FORMAT/AF exists, report it as af. Else, check FORMAT/AD. If not, check FORMAT/*U
    try:
        af = rec.format('AF')[:,0]
    except:
        try:
            alt_counts = rec.format('AD')[:,1]  # AD=REF,ALT so 1 is the position of ALT
        except:
            if rec.is_snp:
                alt_counts = rec.format(rec.ALT[0] + 'U')[:,0]
            else:
                alt_counts = rec.format('TIR')[:,0]
        af = alt_counts / dp

    try:
        mq = rec.format('MQ', float)[:,0]  # VarDict has an incorrect MQ header (with Integer type instead of Float), so need to specify "float" type here explicitly otherwise MQ won't be parsed
    except:
        try:
            mq = rec.format('MMQ')[:,0]
        except:
            mq = [None for _ in [control_index, tumor_index]]
            mq[tumor_index] = rec.INFO.get('MQ')

    return af, dp, mq


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
