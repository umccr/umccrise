from os.path import isfile, join, dirname, abspath
from ngs_utils.file_utils import verify_file
from ngs_utils.key_genes_utils import get_genes_from_file


def package_path():
    return dirname(abspath(__file__))


def get_sig_rmd_file():
    """ Returns path to sig.Rmd file - R-markdown source for mutational signature analysys.
        The file must be located at the same directory as the Snakefile and the patient_analysis module.
    """
    return verify_file(join(package_path(), 'sig.Rmd'))

def get_signatures_probabilities():
    return verify_file(join(package_path(), 'rmd_files', 'signatures_probabilities.txt'), is_critical=True)

def get_suppressors():
    return verify_file(join(package_path(), 'rmd_files', 'suppressors.txt'), is_critical=True)

def get_cancer_genes_ensg():
    return verify_file(join(package_path(), 'ref_data', 'predisposition_genes_engs.txt'), is_critical=True)

def get_key_genes_bed(genome):
    return verify_file(join(package_path(), 'ref_data', 'generated', 'key_genes.' + genome + '.bed'), is_critical=True)

def get_key_genes_set():
    return \
        get_genes_from_file(join(package_path(), 'ref_data', 'az_key_genes.300.txt')) | \
        get_genes_from_file(join(package_path(), 'ref_data', 'umccr_extra_key_genes.txt'))

