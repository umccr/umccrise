#!/usr/bin/env bash

### Cleaning up the environment
unset PYTHONPATH
unset CONDA_PREFIX

GIT_DIR=$(basename $(dirname $(readlink -e $0)))
INSTALL_BASE_DIR=$(realpath $(pwd -P)/miniconda)

### Install conda
if [[ "$OSTYPE" == "darwin"* ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh \
    -O miniconda.sh && chmod +x miniconda.sh
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O miniconda.sh && chmod +x miniconda.sh
fi
bash miniconda.sh -b -p ${INSTALL_BASE_DIR} && rm miniconda.sh
export PATH=${INSTALL_BASE_DIR}/bin:$PATH
conda config --set always_yes yes --set changeps1 no
conda install -c conda-forge 'mamba==0.15.3'
mamba install 'conda==4.11.0'

### Install environments
ENV_NAME=umccrise
mamba env create -p ${INSTALL_BASE_DIR}/envs/${ENV_NAME} --file ${GIT_DIR}/envs/umccrise.yml
mamba env create -p ${INSTALL_BASE_DIR}/envs/${ENV_NAME}_hmf --file ${GIT_DIR}/envs/hmf.yml
mamba env create -p ${INSTALL_BASE_DIR}/envs/${ENV_NAME}_cancer_report --file ${GIT_DIR}/envs/cancer_report.yml
mamba env create -p ${INSTALL_BASE_DIR}/envs/${ENV_NAME}_gatk4 --file ${GIT_DIR}/envs/gatk4.yml
mamba env create -p ${INSTALL_BASE_DIR}/envs/${ENV_NAME}_neoantigens --file ${GIT_DIR}/envs/neoantigens.yml
mamba env create -p ${INSTALL_BASE_DIR}/envs/${ENV_NAME}_oviraptor --file ${GIT_DIR}/envs/oviraptor.yml
if [[ "$OSTYPE" == "darwin"* ]]; then
    mamba env create -p ${INSTALL_BASE_DIR}/envs/${ENV_NAME}_pcgr --file ${GIT_DIR}/envs/pcgr_macos.yml
else
    mamba env create -p ${INSTALL_BASE_DIR}/envs/${ENV_NAME}_pcgr --file ${GIT_DIR}/envs/pcgr_linux.yml
fi

# Instead of `conda activate ${INSTALL_BASE_DIR}/envs/${ENV_NAME}`:
ENV_NAME=umccrise
export PATH=${INSTALL_BASE_DIR}/envs/${ENV_NAME}/bin:$PATH
export CONDA_PREFIX=${INSTALL_BASE_DIR}/envs/${GIT_DIR}

# Install master on top
pip install -e ${GIT_DIR}

### Create the loader script
ENV_NAME=umccrise
cat <<EOT > load_umccrise.sh
unset PYTHONPATH
unset PERL5LIB
export PATH=$PATH
export CONDA_PREFIX=${INSTALL_BASE_DIR}/envs/${ENV_NAME}
EOT

# Clean up
conda clean --yes --tarballs

# TODO: clone vladsaveliev/pVACtools and do `pip install /g/data/gx8/extras/umccrise_2020_Sep/pVACtools`
