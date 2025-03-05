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

RUN mkdir /installed

WORKDIR /installed

## install bowtie2-2.5.4
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip && unzip bowtie2-2.5.4-linux-x86_64.zip && chmod +x bowtie2-2.5.4-linux-x86_64/bowtie2*
ENV PATH=$PATH:/installed/bowtie2-2.5.4-linux-x86_64/

## install humann, metaphlan, and awscli
RUN pip install humann --no-binary :all: 
RUN pip install metaphlan 
RUN pip install awscli

## install SAMtools-1.21
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && tar -xjf samtools-1.21.tar.bz2 && cd samtools-1.21 && make && make install && chmod +x samtools
ENV PATH=$PATH:/installed/samtools-1.21

## install seqtk and fastqc
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
