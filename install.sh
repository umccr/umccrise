#!/usr/bin/env bash

### Clone the repo
git clone https://github.com/umccr/umccrise

### Cleaning up the environment
unset PYTHONPATH
unset CONDA_PREFIX

### Install conda
if [[ "$OSTYPE" == "darwin"* ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh && chmod +x miniconda.sh
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && chmod +x miniconda.sh
fi
./miniconda.sh -b -p $PWD/miniconda && rm miniconda.sh
export PATH=$PWD/miniconda/bin:$PATH
conda update conda

### Install environments
ENV_NAME=umccrise
conda env create -p $PWD/miniconda/envs/${ENV_NAME} --file umccrise/envs/umccrise.yml
conda env create -p $PWD/miniconda/envs/${ENV_NAME}_hmf --file umccrise/envs/hmf.yml
conda env create -p $PWD/miniconda/envs/${ENV_NAME}_multiqc --file umccrise/envs/multiqc.yml
conda env create -p $PWD/miniconda/envs/${ENV_NAME}_cancer_report --file umccrise/envs/cancer_report.yml
if [[ "$OSTYPE" == "darwin"* ]]; then
    conda env create -p $PWD/miniconda/envs/${ENV_NAME}_pcgr --file umccrise/envs/pcgr_macos.yml
else
    conda env create -p $PWD/miniconda/envs/${ENV_NAME}_pcgr --file umccrise/envs/pcgr_linux.yml
fi

# Instead of `conda activate $PWD/miniconda/envs/${ENV_NAME}`:
export PATH=$PWD/miniconda/envs/${ENV_NAME}/bin:$PATH
export CONDA_PREFIX=$PWD/miniconda/envs/umccrise

# Install master on top
pip install -e umccrise

R -e "install.packages('stringi', dependencies = F, repos = 'http://cran.us.r-project.org')"
R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')"
R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')"

### Create the loader script
ENV_NAME=umccrise
cat <<EOT > load_umccrise.sh
unset PYTHONPATH
unset PERL5LIB
export PATH=$PWD/miniconda/envs/${ENV_NAME}/bin:$PWD/miniconda/bin:\$PATH
export CONDA_PREFIX=$PWD/miniconda/envs/${ENV_NAME}
EOT

conda clean --yes --tarballs