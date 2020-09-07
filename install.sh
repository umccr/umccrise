#!/usr/bin/env bash

### Cleaning up the environment
unset PYTHONPATH
unset CONDA_PREFIX

GIT_DIR=$(basename $(dirname $(readlink -e $0)))

### Install conda
if [[ "$OSTYPE" == "darwin"* ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh \
    -O miniconda.sh && chmod +x miniconda.sh
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O miniconda.sh && chmod +x miniconda.sh
fi
./miniconda.sh -b -p $PWD/miniconda && rm miniconda.sh
export PATH=$PWD/miniconda/bin:$PATH
conda config --set always_yes yes --set changeps1 no
conda install -c conda-forge mamba
mamba update conda

### Install environments
ENV_NAME=umccrise
mamba env create -p $PWD/miniconda/envs/${ENV_NAME} --file ${GIT_DIR}/envs/umccrise.yml
mamba env create -p $PWD/miniconda/envs/${ENV_NAME}_hmf --file ${GIT_DIR}/envs/hmf.yml
mamba env create -p $PWD/miniconda/envs/${ENV_NAME}_cancer_report --file ${GIT_DIR}/envs/cancer_report.yml
mamba env create -p $PWD/miniconda/envs/${ENV_NAME}_microbiome --file ${GIT_DIR}/envs/microbiome.yml
if [[ "$OSTYPE" == "darwin"* ]]; then
    mamba env create -p $PWD/miniconda/envs/${ENV_NAME}_pcgr --file ${GIT_DIR}/envs/pcgr_macos.yml
else
    mamba env create -p $PWD/miniconda/envs/${ENV_NAME}_pcgr --file ${GIT_DIR}/envs/pcgr_linux.yml
fi

# Instead of `conda activate $PWD/miniconda/envs/${ENV_NAME}`:
ENV_NAME=umccrise
export PATH=$PWD/miniconda/envs/${ENV_NAME}/bin:$PATH
export CONDA_PREFIX=$PWD/miniconda/envs/${GIT_DIR}

# Install master on top
pip install -e ${GIT_DIR}

### Create the loader script
ENV_NAME=umccrise
cat <<EOT > load_umccrise.sh
unset PYTHONPATH
unset PERL5LIB
export PATH=$PATH
export CONDA_PREFIX=${PWD}/miniconda/envs/${ENV_NAME}
EOT

# Install some missing packages for the cancer report environment
export PATH=${PWD}/miniconda/envs/${ENV_NAME}_cancer_report/bin:$PATH
R -e "install.packages('stringi', dependencies = F, repos = 'http://cran.rstudio.com')"
R -e "install.packages('BiocManager', repos = 'http://cran.rstudio.com')"
R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')"
R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')"

# Clean up
conda clean --yes --tarballs