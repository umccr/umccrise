FROM ubuntu:16.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

# Create non-root user with group
ARG USER=umccrise
ARG GROUP=umccrise
RUN groupadd -g 1000 "${GROUP}" && \
    useradd -g "${GROUP}" -r -m -u 1000 "${USER}" && \
    chmod -R a+rw /home/umccrise/

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget git unzip tar gzip bzip2 g++ make zlib1g-dev nano

# Install fonts for pandoc/rmarkdown
RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections
RUN apt-get install -y ttf-mscorefonts-installer

# Setting locales and timezones, based on https://github.com/jacksoncage/node-docker/blob/master/Dockerfile
# (setting UTC for readr expecting UTC https://rdrr.io/github/tidyverse/readr/src/R/locale.R)
ENV LANGUAGE=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
RUN apt-get update && \
    apt-get install -y locales language-pack-en && \
    locale-gen en_US.UTF-8 && \
    dpkg-reconfigure locales && \
    apt-get install -y tzdata && \
    echo "Etc/UTC" > /etc/timezone && \
    dpkg-reconfigure -f noninteractive tzdata

# Install Umccrise
ADD . umccrise
RUN rm -rf umccrise/.git
RUN bash -xe umccrise/install.sh && \
    chown "${USER}:${GROUP}" load_umccrise.sh && \
    mv load_umccrise.sh "/home/${USER}/"

# Clean up
RUN rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/tmp/* && \
    cd /usr/local && \
    apt-get clean && \
    rm -rf /.cpanm

# Set user and work directory
USER "${USER}"
WORKDIR "/home/${USER}"

# Set environment variables
ENV ENV_NAME=umccrise
ENV HOME="/home/${USER}"
ENV PATH="/miniconda/envs/${ENV_NAME}/bin:${HOME}/miniconda/bin:/home/${USER}/.local/bin:${PATH}"
ENV CONDA_PREFIX="/miniconda/envs/${ENV_NAME}"
