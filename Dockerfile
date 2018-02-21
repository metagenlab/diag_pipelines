FROM cashalow/snakemake

RUN git clone https://github.com/metagenlab/diag_pipelines $pipeline_folder

RUN mkdir /opt/conda/envs/

ENV conda_folder=/opt/conda/envs/

WORKDIR ${main}

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Staphylococcus aureus" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Mycobacterium tuberculosis" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Listeria monocytogenes" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Escherichia coli" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Klebsiella pneumoniae" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Enterococcus faecium" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Acinetobacter baumannii" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Legionella pneumophila" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_enterobase.rules --use-conda --conda-prefix ${conda_folder} --config absolute_path_of_pipeline_folder=/home/pipeline_user/snakemake_pipeline/ logging_folder=/home/pipeline_user/logging/ species="Salmonella enterica" -f all

RUN mkdir -p ${main}/validation/Listeria_monocytogenes/

WORKDIR ${main}/validation/Listeria_monocytogenes/

RUN ln -s /home/pipeline_user/core_genomes/Listeria_monocytogenes/ core_genome

RUN cp /home/pipeline_user/snakemake_pipeline/data/validation_datasets/L_monocytogenes_ASM_NGS_2015/* .

RUN mkdir links
       
RUN snakemake --snakefile targets.rules --use-conda --conda-prefix ${conda_folder} --configfile config.yaml --config ref_ids_for_mapping="" -j 4 

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
