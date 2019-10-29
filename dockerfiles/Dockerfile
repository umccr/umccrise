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

# Install
RUN wget https://raw.githubusercontent.com/umccr/umccrise/master/install.sh
RUN bash -x install.sh

# Clean up
RUN rm -rf umccrise/.git && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/tmp/* && \
    cd /usr/local && \
    apt-get clean && \
    rm -rf /.cpanm
