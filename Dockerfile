FROM continuumio/miniconda3:4.3.27

RUN conda config --add channels defaults

RUN conda config --add channels conda-forge

RUN conda config --add channels bioconda

RUN useradd -r -u 1080 pipeline_user

RUN conda install snakemake=4.8.0=py36_0 fontconfig unzip

ENV main=/home/pipeline_user/

ENV pipeline_folder=${main}/snakemake_pipeline/

RUN git clone https://github.com/metagenlab/diag_pipelines $pipeline_folder

RUN mkdir /opt/conda/envs/

ENV conda_folder=/opt/conda/envs/

WORKDIR /usr/local/bin

RUN wget -qO- https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz > parsnp.tar.gz

RUN tar xf parsnp.tar.gz

RUN mv Parsnp-Linux64-v1.2/parsnp .

RUN rm -rf Parsnp-Linux64-v1.2/

RUN mkdir -p ${main}/data/links

WORKDIR ${main}/data/

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} core_genomes/cgMLST/Staphylococcus_aureus.bed core_genomes/cgMLST/Mycobacterium_tuberculosis.bed core_genomes/cgMLST/Listeria_monocytogenes.bed core_genomes/cgMLST/Klebsiella_pneumoniae.bed core_genomes/cgMLST/Enterococcus_faecium.bed core_genomes/cgMLST/Acinetobacter_baumannii.bed core_genomes/cgMLST/Legionella_pneumophila.bed 

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_enterobase.rules --use-conda --conda-prefix ${conda_folder} core_genomes/cgMLST/Salmonella_enterica.bed core_genomes/cgMLST/Escherichia_coli.bed

RUN cp ${pipeline_folder}/*.tsv . 

RUN cp ${pipeline_folder}/config.yaml .

RUN echo '' > links/Myco-10_S10_L001_R1.fastq.gz

RUN echo '' > links/Myco-10_S10_L001_R2.fastq.gz

RUN mkdir -p core_genomes/parsnp/Mycobacterium_tuberculosis/

RUN echo '' > core_genomes/parsnp/Mycobacterium_tuberculosis/parsnp.xmfa

RUN wget -qO- https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh > ${main}/data/references/mash_sketch.msh

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml  quality/multiqc/assembly/multiqc_report.html samples/M10/resistance/mykrobe.tsv typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/core_parsnp_538048/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/full_genome/M10_assembled_genome/bwa/distances_in_snp_mst_no_st.svg virulence_summary.xlsx typing/mlst/summary.xlsx samples/M10/resistance/m_tuberculosis_resistance_genes_4_db_mutations/summary.xlsx

RUN patch /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch

RUN rm -rf links/

RUN rm core_genomes/parsnp/Staphylococcus_aureus/parsnp.xmfa

RUN rm config.yaml

RUN rm *.tsv

RUN mkdir ${main}/data/analysis/

WORKDIR ${main}/data/analysis/

RUN chown -R pipeline_user ${main}/

USER pipeline_user
