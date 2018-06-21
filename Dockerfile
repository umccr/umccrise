FROM ubuntu:16.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

ENV HOSTNAME umccrise_docker
ENV TEST_DATA_PATH=/umccrise/tests/umccrise_test_data
ENV BCBIO_GENOMES_PATH=/genomes
ENV PON_PATH=/panel_of_normals

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget git unzip tar gzip bzip2 g++ make zlib1g-dev nano

# Setting locales and timezones, based on https://github.com/jacksoncage/node-docker/blob/master/Dockerfile
# (setting UTC for readr expecting UTC https://rdrr.io/github/tidyverse/readr/src/R/locale.R)
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
RUN apt-get update && \
    apt-get install -y locales language-pack-en && \
    locale-gen en_US.UTF-8 && \
    dpkg-reconfigure locales && \
    apt-get install -y tzdata && \
    echo "Etc/UTC" > /etc/timezone && \
    dpkg-reconfigure -f noninteractive tzdata

# Install conda
RUN wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /miniconda

# Instead of `. /miniconda/etc/profile.d/conda.sh`:
ENV PATH /miniconda/bin:$PATH
ENV CONDA_EXE /miniconda/bin/conda
ENV CONDA_ROOT /miniconda
ENV CONDA_PYTHON_EXE /miniconda/bin/python

# Install conda environment
COPY environment.yml .
RUN hash -r && \
    conda config --set always_yes yes --set changeps1 no && \
    conda update -q conda && \
    conda env create -n umccrise --file environment.yml

# Instead of `conda activate umccrise`
ENV PATH /miniconda/envs/umccrise/bin:$PATH
ENV CONDA_PREFIX /miniconda/envs/umccrise
ENV CONDA_DEFAULT_ENV umccrise

# Install source and Peter's circos library
COPY umccrise umccrise/umccrise
COPY scripts umccrise/scripts
COPY vendor umccrise/vendor
COPY tests/test.py umccrise/tests/test.py
COPY setup.py umccrise/setup.py
RUN pip install -e umccrise && \
    R -e "library(devtools) ; options(unzip = '/usr/bin/unzip') ; devtools::install_github('umccr/rock')"

# Clean up
RUN rm -rf umccrise/.git && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/tmp/* && \
    conda clean --yes --tarballs && \
    cd /usr/local && \
    apt-get clean && \
    rm -rf /.cpanm
