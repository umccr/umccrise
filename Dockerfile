FROM ubuntu:21.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

# Create a non-root user and group
ARG USER=umccrise
ARG GROUP=umccrise
RUN groupadd -g 1000 "${GROUP}" && \
    useradd -g "${GROUP}" -r -m -u 1000 "${USER}" && \
    chmod -R a+rw /home/umccrise/

# Retrieve package list, set debconf to non-interactive mode
RUN apt-get update
ARG DEBIAN_FRONTEND='noninteractive'

# Set locale and timezone data
ARG TZ='Etc/UTC'
ENV LANGUAGE='en_US.UTF-8'
ENV LANG='en_US.UTF-8'
ENV LC_ALL='en_US.UTF-8'
RUN apt-get install -y \
        tzdata \
        locales && \
    locale-gen en_US.UTF-8

# Install required base packages
RUN apt-get install -y \
	      git \
	      unzip \
	      wget \
        curl

# Accept the Microsoft TrueType core fonts license and then install non-interactively
# Required for pandas/rmarkdown
RUN echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | \
    debconf-set-selections && \
    apt-get install -y ttf-mscorefonts-installer

# Install Umccrise, remove unused load script
ADD . umccrise
RUN rm -rf umccrise/.git
RUN /bin/bash -xe umccrise/install.sh && \
    rm load_umccrise.sh

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
ENV PATH="/miniconda/envs/${ENV_NAME}/bin:/miniconda/bin:/home/${USER}/.local/bin:${PATH}"
ENV CONDA_PREFIX="/miniconda/envs/${ENV_NAME}"
