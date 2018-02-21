FROM cashalow/snakemake

RUN git clone https://github.com/metagenlab/diag_pipelines $pipeline_folder

RUN mkdir /opt/conda/envs/

ENV conda_folder=/opt/conda/envs/

WORKDIR ${main}

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Staphylococcus aureus"

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Mycobacterium tuberculosis"

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Listeria monocytogenes"

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Escherichia coli"

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Klebsiella pneumoniae"

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Enterococcus faecium"

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Acinetobacter baumannii"

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Legionella pneumophila"

WORKDIR ${main}/data/

RUN mkdir -p ${main}/data/links

RUN cp ${pipeline_folder}/*.tsv . 

RUN cp ${pipeline_folder}/config.yaml .


RUN patch /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch


RUN rm links/*

RUN rm config.yaml

RUN rm *.tsv

RUN rm core_genome/parsnp/parsnp.xmfa

RUN chown pipeline_user -R ${main}

USER pipeline_user
