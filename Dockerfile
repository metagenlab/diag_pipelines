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

COPY ./* /diag_pipeline/
#RUN git clone --single-branch --branch d18793ecba244a3e96a49e511a17106b2e9bd66f https://github.com/metagenlab/diag_pipelines /diag_pipeline

RUN conda install snakemake=5.7.0 singularity=3.0.1

#RUN conda install unzip
#RUN conda install sra-tools=2.9.1
#ENV NCBI_API_KEY=719f6e482d4cdfa315f8d525843c02659408

RUN conda init bash
ENTRYPOINT ["/bin/bash"]
WORKDIR /diag_pipeline/