FROM ubuntu:16.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

ENV HOSTNAME umccrise_docker
ENV TEST_DATA_PATH=/umccrise/umccrise_test_data

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget git unzip tar gzip bzip2 g++ make zlib1g-dev nano

# Install fonts for pandoc/rmarkdown
RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections
RUN apt-get install -y ttf-mscorefonts-installer

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
ENV PATH /miniconda/bin:$PATH

# Install conda environments
COPY envs envs
RUN hash -r && \
    conda config --set always_yes yes --set changeps1 no
RUN conda env create -n umccrise --file envs/umccrise.yml
RUN conda env create -n umccrise_purple --file envs/purple.yml
RUN conda env create -n umccrise_pcgr --file envs/pcgr_linux.yml
# Dirty hack to fix circos https://github.com/bioconda/bioconda-recipes/issues/9830#issuecomment-441438177
RUN cd /miniconda/envs/umccrise_purple/lib && if [ ! -e libwebp.so.7 ] ; then ln -s libwebp.so.6 libwebp.so.7; fi

# Instead of `conda activate umccrise`
ENV PATH /miniconda/envs/umccrise/bin:$PATH
ENV CONDA_PREFIX /miniconda/envs/umccrise

# Clean up
RUN rm -rf umccrise/.git && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/tmp/* && \
    conda clean --yes --tarballs && \
    cd /usr/local && \
    apt-get clean && \
    rm -rf /.cpanm

# Install source
COPY umccrise umccrise/umccrise
COPY scripts umccrise/scripts
COPY vendor umccrise/vendor
COPY setup.py umccrise/setup.py
RUN pip install -e umccrise

RUN R -e "install.packages('stringi', dependencies = F, repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')"
RUN R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')"
