FROM continuumio/miniconda3

RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN useradd -r -u 1080 pipeline_user

RUN apt-get install -y fontconfig unzip

RUN conda install snakemake=4.8.0=py36_0

ENV main=/home/pipeline_user/

ENV pipeline_folder=${main}/snakemake_pipeline/

RUN git clone https://github.com/metagenlab/diag_pipelines $pipeline_folder

RUN mkdir /opt/conda/envs/

ENV conda_folder=/opt/conda/envs/

RUN mkdir -p ${main}/data/links

WORKDIR ${main}/data/

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config species="Staphylococcus aureus" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config species="Mycobacterium tuberculosis" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config species="Listeria monocytogenes" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config species="Escherichia coli" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config species="Klebsiella pneumoniae" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config species="Enterococcus faecium" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config species="Acinetobacter baumannii" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} --config species="Legionella pneumophila" -f all

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_enterobase.rules --use-conda --conda-prefix ${conda_folder} --config species="Salmonella enterica" -f all

RUN cp ${pipeline_folder}/*.tsv . 

RUN cp ${pipeline_folder}/config.yaml .

RUN echo '' > links/Staaur-10_S10_L001_R1.fastq.gz

RUN echo '' > links/Staaur-10_S10_L001_R2.fastq.gz

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml  quality/multiqc/assembly/multiqc_report.html samples/S10/resistance/mykrobe.tsv samples/S10/resistance/rgi.tsv typing/freebayes_joint_genotyping/core_ridom/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/core_ridom/bwa/distances_in_snp_mst_no_st.svg

RUN patch /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch

RUN rm -rf links/

RUN rm config.yaml

RUN rm *.tsv

WORKDIR /usr/local/bin

RUN wget -qO- https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz > parsnp.tar.gz

RUN tar xf parsnp.tar.gz

RUN mv Parsnp-Linux64-v1.2/parsnp .

RUN rm -rf Parsnp-Linux64-v1.2/

WORKDIR ${main}/data/

RUN chown -R pipeline_user ${main}/

USER pipeline_user



