FROM ubuntu:16.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

ENV HOSTNAME umccrise_docker
ENV TEST_DATA_PATH=/umccrise/tests/umccrise_test_data
ENV BCBIO_GENOMES_PATH=/genomes
ENV PON_PATH=/panel_of_normals

VOLUME $TEST_DATA_PATH
VOLUME $BCBIO_GENOMES_PATH
VOLUME $PON_PATH

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

#RUN apt-get update && \
#    apt-get install -y make g++
#    apt-get upgrade -y libstdc++6
#-t vieux/apache:2.0

# Install environemnt
COPY environment.yml .
RUN conda env create -n umccrise --file environment.yml
RUN conda info -a

# Instead of `conda activate umccrise`:
ENV PATH /miniconda/envs/umccrise/bin:$PATH
ENV CONDA_PREFIX /miniconda/envs/umccrise
ENV CONDA_DEFAULT_ENV umccrise

# Copy and install source
COPY umccrise umccrise/umccrise
COPY scripts umccrise/scripts
COPY vendor umccrise/vendor
COPY tests/test.py umccrise/tests/test.py
COPY setup.py umccrise/setup.py
RUN pip install -e umccrise
COPY /Users/vsaveliev/git/umccr/python_utils python_utils
RUN pip install -e python_utils

# Clean up
RUN rm -rf umccrise/.git && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/tmp/* && \
    conda clean --yes --tarballs && \
    cd /usr/local && \
    apt-get clean && \
    rm -rf /.cpanm

