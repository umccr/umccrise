#!/usr/bin/env bash

### Clone the repo
git clone https://github.com/umccr/umccrise

### Install conda
# source ~/reset_path.sh  # try to clean up your PATH variable to avoid messing with other environments
unset PYTHONPATH
unset CONDA_PREFIX
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda && rm miniconda.sh
export PATH=$(pwd)/miniconda/bin:$PATH
conda update conda

### Install environments
ENV_NAME=umccrise
conda env create -p $(pwd)/miniconda/envs/${ENV_NAME} --file umccrise/envs/umccrise.yml
conda env create -p $(pwd)/miniconda/envs/${ENV_NAME}_purple --file umccrise/envs/purple.yml
# conda env create -p $(pwd)/miniconda/envs/${ENV_NAME}_pcgr --file umccrise/envs/pcgr_macos.yml    # macos
conda env create -p $(pwd)/miniconda/envs/${ENV_NAME}_pcgr --file umccrise/envs/pcgr_linux.yml    # linux
conda activate $(pwd)/miniconda/envs/${ENV_NAME}  # export PATH=$(pwd)/miniconda/envs/${ENV_NAME}/bin:$PATH
pip install -e umccrise

### Create the loader script
ENV_NAME=umccrise
cat <<EOT > load_umccrise.sh
unset PYTHONPATH
unset PERL5LIB
export PATH=$(pwd)/miniconda/envs/${ENV_NAME}/bin:$(pwd)/miniconda/bin:\$PATH
export CONDA_PREFIX=$(pwd)/miniconda/envs/${ENV_NAME}
EOT

