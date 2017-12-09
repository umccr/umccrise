from os.path import isfile, join, dirname, abspath


def get_sig_rmd_file():
    """ Returns path to sig.Rmd file - R-markdown source for mutational signature analysys.
        The file must be located at the same directory as the Snakefile and the patient_analysis module.
    """
    sig_rmd = join(dirname(abspath(__file__)), 'sig.Rmd')
    assert isfile(sig_rmd)
    return sig_rmd
