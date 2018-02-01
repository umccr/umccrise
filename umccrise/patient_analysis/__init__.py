from os.path import isfile, join, dirname, abspath
from ngs_utils.file_utils import verify_file

here = dirname(abspath(__file__))


def package_path():
    return here


def get_sig_rmd_file():
    """ Returns path to sig.Rmd file - R-markdown source for mutational signature analysys.
        The file must be located at the same directory as the Snakefile and the patient_analysis module.
    """
    return verify_file(join(here, 'sig.Rmd'))

def get_signatures_probabilities():
    return verify_file(join(here, 'rmd_files', 'signatures_probabilities.txt'))

def get_suppressors():
    return verify_file(join(here, 'rmd_files', 'suppressors.txt'))

def get_cancer_genes_ensg():
    return verify_file(join(here, 'cancer_genes_ENSG.txt'))
