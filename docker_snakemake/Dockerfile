FROM continuumio/miniconda3

RUN /bin/bash -c "conda config --add channels defaults"
RUN /bin/bash -c "conda config --add channels conda-forge"
RUN /bin/bash -c "conda config --add channels bioconda"

RUN useradd -r -u 1080 pipeline_user

RUN apt-get install -y fontconfig unzip

RUN conda install snakemake=4.6.0=py36_0

ENV main=/home/pipeline_user/

ENV pipeline_folder=${main}/snakemake_pipeline/

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

RUN mkdir -p ${main}/data/links

WORKDIR ${main}/data/

RUN cp ${pipeline_folder}/*.tsv . 

RUN cp ${pipeline_folder}/config.yaml .

RUN echo '' > links/Staaur-10_S10_L001_R1.fastq.gz

RUN echo '' > links/Staaur-10_S10_L001_R2.fastq.gz

RUN ln -s ${main}/core_genomes/Staphylococcus_aureus/ core_genome

RUN snakemake --snakefile ${pipeline_folder}/workflows/assembly_quality.rules --config logging_folder="/home/pipeline_user/logging/" --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml  quality/multiqc/self_genome/multiqc_report.html 

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --config logging_folder="/home/pipeline_user/logging/" --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml resistance_summary.xlsx

RUN snakemake --snakefile ${pipeline_folder}/workflows/virulence.rules --config logging_folder="/home/pipeline_user/logging/" --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml virulence_summary.xlsx

RUN snakemake --snakefile ${pipeline_folder}/workflows/typing.rules --configfile config.yaml --config species="Staphylococcus aureus" logging_folder="/home/pipeline_user/logging/" --use-conda --create-envs-only --conda-prefix ${conda_folder}  typing/freebayes_joint_genotyping/core_ridom/33148/bwa/distances_snp_mst.pdf

RUN snakemake --snakefile ${pipeline_folder}/workflows/typing.rules --config taxid=1280 logging_folder="/home/pipeline_user/logging/" --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml typing/gatk_gvcfs/core_ridom/33148/bwa/distances_snp_mst.pdf

RUN patch /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch

RUN rm links/*

RUN rm config.yaml

RUN rm *.tsv



