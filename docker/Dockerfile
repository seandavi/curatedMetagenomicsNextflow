FROM continuumio/miniconda3

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y curl
RUN apt-get install -y wget

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

## install metaphlan3, humann3, and the curatedMetagenomics pipeline
RUN conda update conda &&\
    conda install -c conda-forge awscli &&\
    conda install -c bioconda metaphlan=3 &&\
    conda install -c biobakery humann=3 &&\
    conda install fastqc &&\
    conda install seqtk &&\
    conda install -c bioconda fastq-scan

## install sratoolkit
RUN mkdir /tmp/sra
WORKDIR /tmp/sra

RUN curl -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/setup-apt.sh
RUN /bin/bash setup-apt.sh
ENV PATH=/usr/local/ncbi/sra-tools/bin:$PATH

## interactive user settings workaround
COPY usersettings /root/.ncbi/user-settings.mkfg

RUN rm -rf /tmp/sra

## display versions
RUN fastq-dump --version
RUN metaphlan --version
RUN humann --version

# install ps for nextflow task monitoring
RUN apt-get update && apt-get install -y procps pigz pv

## root directory for databases
ENV mpa_dir=/usr/local/miniconda3/lib/python3.7/site-packages/metaphlan
ENV hnn_dir=/usr/local/miniconda3/lib/python3.7/site-packages/humann
WORKDIR /root

CMD ["/bin/bash"]
