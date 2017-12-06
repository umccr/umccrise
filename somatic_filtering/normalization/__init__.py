import os
from ngs_utils.file_utils import verify_file
from os.path import isfile
from somatic_filtering.utils import get_ref_path

def make_normalize_cmd(input_file, output_file, reference_fasta):
    if not isfile(reference_fasta):
        reference_fasta = os.path.join(get_ref_path(), reference_fasta)
    verify_file(reference_fasta, is_critical=True)
    return (
        f'bcftools norm -m \'-\' {input_file} -Ov -f {reference_fasta}'
        f' | vcfallelicprimitives -t DECOMPOSED --keep-geno | vcfstreamsort | bgzip -c > {output_file}'
        f' && tabix -p vcf {output_file}')

