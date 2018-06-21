FROM ubuntu:16.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget git unzip tar gzip bzip2 g++ make zlib1g-dev nano

# Install conda
RUN wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /miniconda

# Instead of `. /miniconda/etc/profile.d/conda.sh`
ENV PATH /miniconda/bin:$PATH
ENV CONDA_EXE /miniconda/bin/conda
ENV CONDA_ROOT /miniconda
ENV CONDA_PYTHON_EXE /miniconda/bin/python

RUN hash -r && \
    conda config --set always_yes yes --set changeps1 no && \
    conda update -q conda

# Install environemnt
COPY environment.yml .
RUN conda env create -n umccrise --file environment.yml
RUN conda info -a

# Instead of `conda activate umccrise`:
ENV PATH /miniconda/envs/umccrise/bin:$PATH
ENV CONDA_PREFIX /miniconda/envs/umccrise
ENV CONDA_DEFAULT_ENV umccrise

# Install Peter's circos library
RUN R -e "library(devtools) ; options(unzip = '/usr/bin/unzip') ; devtools::install_github('umccr/rock')"

# Setting locales and timezones, based on https://github.com/jacksoncage/node-docker/blob/master/Dockerfile
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
RUN apt-get update && \
    apt-get install -y locales && \
    locale-gen en_US.UTF-8 && \
    dpkg-reconfigure locales
# (setting UTC for readr expecting UTC https://rdrr.io/github/tidyverse/readr/src/R/locale.R)
RUN apt-get update && \
    apt-get install -y tzdata && \
    echo "Etc/UTC" > /etc/timezone && \
    dpkg-reconfigure -f noninteractive tzdata

# Copy and install source
COPY umccrise umccrise/umccrise
COPY scripts umccrise/scripts
COPY vendor umccrise/vendor
COPY tests/test.py umccrise/tests/test.py
COPY setup.py umccrise/setup.py
RUN pip install -e umccrise

# Clean up
RUN rm -rf umccrise/.git && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/tmp/* && \
    conda clean --yes --tarballs && \
    cd /usr/local && \
    apt-get clean && \
    rm -rf /.cpanm

ENV HOSTNAME umccrise_docker
ENV TEST_DATA_PATH=/umccrise/tests/umccrise_test_data
ENV BCBIO_GENOMES_PATH=/genomes
ENV PON_PATH=/panel_of_normals
