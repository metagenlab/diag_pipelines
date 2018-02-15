FROM continuumio/miniconda3

RUN /bin/bash -c "conda config --add channels defaults"
RUN /bin/bash -c "conda config --add channels conda-forge"
RUN /bin/bash -c "conda config --add channels bioconda"

RUN useradd -r -u 1080 pipeline_user

RUN apt-get install -y fontconfig

ENV main=/home/pipeline_user/

ENV pipeline_folder=${main}/snakemake_pipeline/

RUN git clone --branch v0.1.3 https://github.com/metagenlab/diag_pipelines $pipeline_folder

RUN conda install snakemake=4.6.0=py36_0

RUN mkdir -p ${main}/data/links

WORKDIR ${main}/data/

RUN cp ${pipeline_folder}/*.tsv . 

RUN cp ${pipeline_folder}/config.yaml .

RUN echo '' > links/Staaur-10_S10_L001_R1.fastq.gz

RUN echo '' > links/Staaur-10_S10_L001_R2.fastq.gz

RUN mkdir -p core_genome/parsnp

RUN echo '' > core_genome/parsnp/parsnp.xmfa

RUN mkdir /opt/conda/envs/

ENV conda_folder=/opt/conda/envs/

RUN snakemake --snakefile ${pipeline_folder}/workflows/ring_trial/target_files.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml

RUN patch /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch

RUN rm links/*

RUN rm config.yaml

RUN rm *.tsv

RUN rm core_genome/parsnp/parsnp.xmfa

RUN chown pipeline_user -R ${main}

USER pipeline_user
