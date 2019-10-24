#!/usr/bin/env bash

### Clone the repo
git clone --recursive https://github.com/umccr/umccrise

### Install conda
# source ~/reset_path.sh  # try to clean up your PATH variable to avoid messing with other environments
unset PYTHONPATH
unset CONDA_PREFIX
if [[ "$OSTYPE" == "darwin"* ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
fi
bash miniconda.sh -b -p ./miniconda && rm miniconda.sh
export PATH=$(pwd)/miniconda/bin:$PATH
conda update conda

### Install environments
ENV_NAME=umccrise
conda env create -y -p $(pwd)/miniconda/envs/${ENV_NAME} --file umccrise/envs/umccrise.yml
conda env create -y -p $(pwd)/miniconda/envs/${ENV_NAME}_hmf --file umccrise/envs/hmf.yml
if [[ "$OSTYPE" == "darwin"* ]]; then
    conda env create -y -p $(pwd)/miniconda/envs/${ENV_NAME}_pcgr --file umccrise/envs/pcgr_macos.yml
else
    conda env create -y -p $(pwd)/miniconda/envs/${ENV_NAME}_pcgr --file umccrise/envs/pcgr_linux.yml
fi

# Instead of `conda activate $(pwd)/miniconda/envs/${ENV_NAME}`:
export PATH=$(pwd)/miniconda/envs/${ENV_NAME}/bin:$PATH
export CONDA_PREFIX=$(pwd)/miniconda/envs/umccrise

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
export PATH=$(pwd)/miniconda/envs/${ENV_NAME}/bin:$(pwd)/miniconda/bin:\$PATH
export CONDA_PREFIX=$(pwd)/miniconda/envs/${ENV_NAME}
EOT
