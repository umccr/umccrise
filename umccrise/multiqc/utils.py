from __future__ import division
import os


def format_decimal(value, unit=None):
    if value is None:
        return None
    unit_str = ('<span style="font-size: 50%; line-height: 1;">&nbsp;</span>' + unit) if unit else ''
    if value <= 9999:
        return '{value}{unit_str}'.format(**locals())
    else:
        v = '{value:,}{unit_str}'.format(**locals())
        return v.replace(',', '<span style="margin-left: 0.2em;"></span>')


def format_bed_info(d, genome_info):
    bed, size, regions, genes = d['bed'], \
                                format_decimal(d.get('size')), \
                                format_decimal(d.get('regions')), \
                                format_decimal(d.get('genes'))
    bed_name = os.path.basename(bed)
    html = '{bed_name}'
    if size is not None:
        percent = (100.0 * d['size'] / genome_info['size']) if genome_info.get('size') else 0
        html += ' ({size} bp'
        if percent >= 0.01:
            html += ' or {percent:.2f}% of the genome'
        if regions is not None:
            html += ', {regions} region' + ('s' if d.get('regions') != 1 else '')
        if genes is not None:
            html += ', {genes} gene' + ('s' if d.get('genes') != 1 else '')
        html += ')'
    return html.format(**locals())
