# Base image for every process without its own container, and the container
# of the default MetaPhlAn profile (conf/metaphlan_profiles.config). Published
# to ghcr.io/seandavi/curatedmetagenomics:metaphlan<version> by
# .github/workflows/base-image.yml. No HUMAnN: bundles bring their own (#89).
#
# bookworm, not the floating python:3.9 tag (now trixie), keeps the OS
# toolchain (Java 17, pigz, pv) the same as the metaphlan4.2.2 image.
FROM python:3.9-bookworm

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y git make curl wget gcc g++ unzip

# install ps for nextflow task monitoring
RUN apt-get update && apt-get install -y procps pigz pv default-jre

## install sratoolkit
RUN mkdir /usr/local/ncbi
WORKDIR /usr/local/ncbi

# Pinned: the PATH below names the version, and sdk/current has moved past it.
RUN curl -o sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz
RUN tar -vxzf sratoolkit.tar.gz
ENV PATH=/usr/local/ncbi/sratoolkit.3.2.1-ubuntu64/bin:$PATH

## interactive user settings workaround
COPY usersettings /root/.ncbi/user-settings.mkfg

RUN rm -rf /tmp/sra

RUN mkdir /installed

WORKDIR /installed

## install bowtie2-2.5.4
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip && unzip bowtie2-2.5.4-linux-x86_64.zip && chmod +x bowtie2-2.5.4-linux-x86_64/bowtie2*
ENV PATH=$PATH:/installed/bowtie2-2.5.4-linux-x86_64/

## install metaphlan
# 4.2.6 is a GitHub release only (not on PyPI); bioconda builds the same tarball.
# Its --version still reports 4.2.5 (upstream did not bump __version__).
RUN pip install "https://github.com/biobakery/MetaPhlAn/archive/refs/tags/4.2.6.tar.gz#sha256=b72991c2d7d47818c99959f51b4bc4137fe0bfbf8b3d76e1c2e895b499d6e54c"

## install awscli
RUN pip install awscli

## install SAMtools-1.21
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && tar -xjf samtools-1.21.tar.bz2 && cd samtools-1.21 && make && make install && chmod +x samtools
ENV PATH=$PATH:/installed/samtools-1.21

## install seqtk and fastqc
RUN git clone https://github.com/lh3/seqtk.git && cd seqtk && make
RUN chmod +x seqtk/seqtk
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && unzip fastqc_v0.12.1.zip && chmod +x FastQC/fastqc
ENV PATH=$PATH:/installed/seqtk:/installed/FastQC/

## install minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30_x64-linux.tar.bz2 | tar -jxvf -
ENV PATH=$PATH:/installed/minimap2-2.30_x64-linux/

## install kneaddata
RUN pip install kneaddata==0.12.3 --install-option='--bypass-dependencies-install'

## install trimmomatic for kneaddata
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && unzip Trimmomatic-0.39.zip

WORKDIR /installed

## display versions
RUN echo 'VERSIONS' > versions.txt
RUN echo 'fastq-dump' >> versions.txt
RUN fastq-dump --version >> versions.txt
RUN echo 'metaphlan' >> versions.txt
RUN metaphlan --version >> versions.txt

## root directory for databases
ENV mpa_dir=/usr/local/miniconda3/lib/python3.7/site-packages/metaphlan
WORKDIR /root

CMD ["/bin/bash"]
