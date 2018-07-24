from os.path import isfile, join, dirname, abspath
from ngs_utils.file_utils import verify_file


def package_path():
    return dirname(abspath(__file__))


def get_sig_rmd_file():
    """ Returns path to sig.Rmd file - R-markdown source for mutational signature analysys.
        The file must be located at the same directory as the Snakefile and the patient_analysis module.
    """
    return verify_file(join(package_path(), 'rmd_files', 'sig.Rmd'), is_critical=True)

def get_signatures_probabilities():
    return verify_file(join(package_path(), 'rmd_files', 'signatures_probabilities.txt'), is_critical=True)
