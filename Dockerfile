FROM ubuntu:16.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget git unzip tar gzip bzip2 g++ make zlib1g-dev

# Install conda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /miniconda

ENV PATH "/miniconda/bin:$PATH"

RUN hash -r && \
    conda config --set always_yes yes --set changeps1 no && \
    conda update -q conda

#RUN apt-get update && \
#    apt-get install -y make g++
#    apt-get upgrade -y libstdc++6
#-t vieux/apache:2.0

RUN conda install -c bioconda -y pybedtools

# Install environemnt
COPY environment.yml .
RUN conda env create -n umccrise --file environment.yml && \
    source activate umccrise && \
    conda info -a

# Download the reference data
RUN wget --no-check-certificate -c https://s3.amazonaws.com/biodata/genomes/GRCh37-seq.tar.gz && \
    tar -xzvpf GRCh37-seq.tar.gz --directory / && \
    gunzip -c /seq/GRCh37.fa.gz > /seq/GRCh37.fa

# Clone the test data
RUN git clone https://github.com/umccr/umccrise_test_data tests/umccrise_test_data && \
    ln -s /seq tests/umccrise_test_data/data/genomes/Hsapiens/GRCh37/seq

# Copy and install source
COPY . umccrise
RUN pip install -e umccrise

RUN apt-get install nano

# Clean up
RUN rm -rf umccrise/.git && \
    rm -rf /var/lib/apt/lists/* /var/tmp/* && \
    conda clean --yes --tarballs && \
    cd /usr/local && apt-get clean && \
    rm -rf /.cpanm

