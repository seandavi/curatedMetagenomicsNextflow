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
RUN conda update conda
RUN conda install -c bioconda metaphlan=3
RUN conda install -c biobakery humann=3
RUN conda install bowtie2=2.3.5.1

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

WORKDIR /usr/local
RUN wget -O bowtie2-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/bowtie2-2.4.1-linux-x86_64.zip/download
RUN apt-get install -y unzip
RUN unzip bowtie2-linux-x86_64.zip
ENV PATH=/usr/local/bowtie2-2.4.1-linux-x86_64:$PATH

# install ps for nextflow task monitoring
RUN apt-get update && apt-get install -y procps

## root directory for databases
ENV mpa_dir=/usr/local/miniconda3/lib/python3.7/site-packages/metaphlan
ENV hnn_dir=/usr/local/miniconda3/lib/python3.7/site-packages/humann
WORKDIR /root

CMD ["/bin/bash"]
