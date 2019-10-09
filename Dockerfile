# Base Image
FROM continuumio/miniconda3:4.7.10
 
################## METADATA ######################
 
LABEL base.image="biocontainers:latest"
LABEL version="1"
LABEL software="diag-pipelines-singularity"
LABEL software.version="1.0"
LABEL description="combined toolset analysis"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Trestan Pillonel

################## INSTALLATION ######################
RUN conda config --add channels defaults  && conda config --add channels conda-forge && conda config --add channels bioconda

RUN mkdir -p /diag_pipeline

ADD ./data /diag_pipeline/data/
ADD ./workflows /diag_pipeline/workflows/
ADD ./rules /diag_pipeline/rules/
ADD ./envs /diag_pipeline/envs/
#RUN git clone --single-branch --branch d18793ecba244a3e96a49e511a17106b2e9bd66f https://github.com/metagenlab/diag_pipelines /diag_pipeline

RUN conda install snakemake=5.7.0

RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools


WORKDIR /usr/local/bin/
RUN wget https://dl.google.com/go/go1.13.1.linux-amd64.tar.gz
RUN tar -C /usr/local -xzf go1.13.1.linux-amd64.tar.gz
ENV GOPATH=/usr/local/bin/
ENV PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin

RUN mkdir -p $GOPATH/src/github.com/sylabs
WORKDIR $GOPATH/src/github.com/sylabs
RUN git clone https://@github.com/sylabs/singularity.git
RUN mkdir -p /usr/local/bin/bin
WORKDIR $GOPATH/src/github.com/sylabs/singularity
RUN apt-get install curl -y
RUN curl https://raw.githubusercontent.com/golang/dep/master/install.sh | sh
RUN ./mconfig
RUN make -C builddir
RUN make -C builddir install

#RUN conda install unzip
#RUN conda install sra-tools=2.9.1
#ENV NCBI_API_KEY=719f6e482d4cdfa315f8d525843c02659408

RUN conda init bash
ENTRYPOINT ["/bin/bash"]
WORKDIR /diag_pipeline/