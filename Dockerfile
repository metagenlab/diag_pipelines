FROM continuumio/miniconda3

RUN /bin/bash -c "conda config --add channels defaults"
RUN /bin/bash -c "conda config --add channels conda-forge"
RUN /bin/bash -c "conda config --add channels bioconda"

RUN useradd -r -u 1080 pipeline_user

RUN apt-get install -y fontconfig

RUN git clone https://github.com/metagenlab/diag_pipelines /snakemake_pipeline/

RUN conda install snakemake=4.6.0=py36_0

RUN mkdir /opt/conda/envs/

ENV conda_folder=/opt/conda/envs/

ENV pipeline_folder=/snakemake_pipeline/

ENV main=/home/pipeline_user/

WORKDIR ${main}

RUN mkdir -p links/

RUN cp ${pipeline_folder}/*.tsv . 

RUN cp ${pipeline_folder}/config.yaml .

RUN echo '' > links/Staaur-10_S10_L001_R1.fastq.gz

RUN echo '' > links/Staaur-10_S10_L001_R2.fastq.gz

RUN mkdir -p core_genome/parsnp

RUN echo '' > core_genome/parsnp/parsnp.xmfa

RUN snakemake --snakefile ${pipeline_folder}/workflows/ring_trial/target_files.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml

RUN patch /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch

RUN rm links/*

RUN rm config.yaml

RUN rm *.tsv

RUN rm core_genome/parsnp/parsnp.xmfa

USER pipeline_user
