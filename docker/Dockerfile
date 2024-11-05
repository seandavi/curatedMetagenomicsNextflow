FROM --platform=linux/amd64 python:3.9

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y git make curl wget gcc g++ unzip


## install sratoolkit
RUN mkdir /tmp/sra
WORKDIR /tmp/sra

RUN curl -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/setup-apt.sh
RUN /bin/bash setup-apt.sh
ENV PATH=/usr/local/ncbi/sra-tools/bin:$PATH

## interactive user settings workaround
COPY usersettings /root/.ncbi/user-settings.mkfg

RUN rm -rf /tmp/sra

RUN pip install humann --no-binary :all: 
RUN pip install metaphlan 
RUN pip install awscli

RUN mkdir /installed

WORKDIR /installed

## install seqtk
RUN git clone https://github.com/lh3/seqtk.git && cd seqtk && make
RUN chmod +x seqtk/seqtk
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && unzip fastqc_v0.12.1.zip && chmod +x FastQC/fastqc
ENV PATH=$PATH:/installed/seqtk:/installed/FastQC/

## install kneaddata
RUN pip install kneaddata --install-option='--bypass-dependencies-install'

## install trimmomatic for kneaddata
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && unzip Trimmomatic-0.39.zip

## display versions
RUN fastq-dump --version
RUN metaphlan --version
RUN humann --version

# install ps for nextflow task monitoring
RUN apt-get update && apt-get install -y procps pigz pv default-jre

## root directory for databases
ENV mpa_dir=/usr/local/miniconda3/lib/python3.7/site-packages/metaphlan
ENV hnn_dir=/usr/local/miniconda3/lib/python3.7/site-packages/humann
WORKDIR /root

CMD ["/bin/bash"]
