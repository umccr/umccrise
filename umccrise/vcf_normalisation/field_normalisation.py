import sys

from ngs_utils.file_utils import add_suffix, get_ungz_gz
from cyvcf2 import VCF, Writer


def use_pysam(vcf_file, proc_fields):
    """ Working.
        Time:           3:34.04
    """
    import pysam
    vcf = pysam.VariantFile(vcf_file)
    vcf.header.filters.add('MSI_FAIL', None, None, 'Possible homopolymer artefact')
    sys.stdout.write(str(vcf.header))
    for rec in vcf:
        msi_fail = proc_fields(rec.ref, rec.alts[0], rec.samples.values()[0]['AF'], rec.info['MSI'])
        if msi_fail:
            rec.filter.add('MSI_FAIL')
        sys.stdout.write(str(rec))

def use_cyvcf2(vcf_file, proc_fields, vcf_out=None):
    """ Working.
        File out:       2:17.51
        stdout + bgzip: 2:50.35
    """
    vcf = VCF(vcf_file)
    vcf.add_filter_to_header({'ID': 'MSI_FAIL', 'Description': 'Possible homopolymer artefact'})
    if vcf_out:
        w = Writer(vcf_out, vcf)
    else:
        w = None
        sys.stdout.write(vcf.raw_header)
    for rec in vcf:
        msi_fail = proc_fields(rec.REF, rec.ALT[0], rec.format('AF')[0][0], rec.INFO['MSI'])
        if msi_fail:
            filters = rec.FILTER.split(';') if rec.FILTER else []
            filters.append('MSI_FAIL')
            rec.FILTER = ';'.join(filters)
        if w:
            w.write_record(rec)
        else:
            sys.stdout.write(str(rec))
    if w:
        w.close()

def main():
    vcf_file = sys.argv[1]
    cmd = sys.argv[2]

    if cmd == 'pysam':
        use_pysam(vcf_file)
    if cmd == 'cyvcf2':
        vcf_out = sys.argv[3] if len(sys.argv) > 3 else None
        use_cyvcf2(vcf_file, vcf_out)


if __name__ == '__main__':
    main()
