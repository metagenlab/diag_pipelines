FROM continuumio/miniconda3

USER pipeline

RUN /bin/bash -c "conda config --add channels defaults"
RUN /bin/bash -c "conda config --add channels conda-forge"
RUN /bin/bash -c "conda config --add channels bioconda"

RUN apt-get install fontconfig

RUN mkdir pipeline/

RUN git clone https://github.com/metagenlab/diag_pipelines /home/pipeline/

RUN conda install snakemake=4.6.0=py36_0

RUN mkdir /opt/conda/envs/

ENV conda_folder=/opt/conda/envs/

ENV pipeline_folder=/home/pipeline/

WORKDIR /usr/local/bin

RUN wget https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz

RUN tar -xf parsnp-Linux64-v1.2.tar.gz

RUN mv Parsnp-Linux64-v1.2/parsnp parsnp

RUN rm -rf Parsnp-Linux64-v1.2
	
RUN rm parsnp-Linux64-v1.2.tar.gz*

RUN mkdir -p /data/links/

WORKDIR /data/ 

RUN cp ${pipeline_folder}/*.tsv . 

RUN cp ${pipeline_folder}/config.yaml .

RUN echo '' > example_sras.tsv

RUN echo '' > links/Staaur-10_S10_L001_R1.fastq.gz

RUN echo '' > links/Staaur-10_S10_L001_R2.fastq.gz

RUN snakemake --snakefile ${pipeline_folder}/workflows/general_workflow.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml

RUN patch /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch

RUN rm links/*
