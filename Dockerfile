FROM continuumio/miniconda3:4.3.27

RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

RUN useradd -r -u 1080 pipeline_user

RUN apt-get install -y fontconfig

RUN conda install snakemake=5.3.0

RUN conda install unzip

RUN conda install sra-tools=2.9.0

#RUN conda install biopython=1.70

ENV main=/home/pipeline_user

ENV pipeline_folder=${main}/snakemake_pipeline

RUN git clone https://github.com/metagenlab/diag_pipelines $pipeline_folder

RUN mkdir /opt/conda/envs/

ENV conda_folder=/opt/conda/envs/

WORKDIR /usr/local/bin

RUN wget -qO- https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz > parsnp.tar.gz

RUN tar xf parsnp.tar.gz && mv Parsnp-Linux64-v1.2/parsnp . && rm -rf Parsnp-Linux64-v1.2/

RUN mkdir -p ${main}/data/links

WORKDIR ${main}/data/

ENV NCBI_API_KEY=719f6e482d4cdfa315f8d525843c02659408

RUN vdb-config --restore-defaults

# setup VFDB database

RUN snakemake --snakefile ${pipeline_folder}/rules/downloading/fetch_VFDB.rules --use-conda --conda-prefix ${conda_folder} references/virulence/VFDB_larger_50aa.faa.phr references/virulence/VFDB_annotations.tab -j 4

#CREATE cgMLST BED FILES

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} core_genomes/cgMLST/Staphylococcus_aureus.bed core_genomes/cgMLST/Mycobacterium_tuberculosis.bed core_genomes/cgMLST/Listeria_monocytogenes.bed core_genomes/cgMLST/Klebsiella_pneumoniae.bed core_genomes/cgMLST/Enterococcus_faecium.bed core_genomes/cgMLST/Acinetobacter_baumannii.bed core_genomes/cgMLST/Legionella_pneumophila.bed core_genomes/cgMLST/Clostridioides_difficile.bed -j 4

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_enterobase.rules --use-conda --conda-prefix ${conda_folder} core_genomes/cgMLST/Salmonella_enterica.bed core_genomes/cgMLST/Escherichia_coli.bed references/538048/genome_gbwithparts.gbwithparts -j 4

RUN snakemake --snakefile ${pipeline_folder}/rules/downloading/adapting_genome_files.rules --use-conda --conda-prefix ${conda_folder} references/538048/genome_gbwithparts.gbk

# CREATE MTB RESISTANCE DATABASES

RUN mkdir -p resistance_db/Mycobacterium_tuberculosis/mutations/

RUN mkdir -p resistance_db/Mycobacterium_tuberculosis/metadata/

RUN cp ${pipeline_folder}/data/Mycobacterium_tuberculosis/mutations/* resistance_db/Mycobacterium_tuberculosis/mutations/

RUN cp ${pipeline_folder}/data/Mycobacterium_tuberculosis/metadata/* resistance_db/Mycobacterium_tuberculosis/metadata/

RUN snakemake --snakefile ${pipeline_folder}/workflows/check_resistance_databases.rules --use-conda --conda-prefix ${conda_folder} resistance_db/Mycobacterium_tuberculosis/mutations/currated_db_all/correct.bed resistance_db/Mycobacterium_tuberculosis/mutations/level_four_agreement/correct.bed resistance_db/Mycobacterium_tuberculosis/mutations/rgi_annotated_full_2_0_0/correct.bed resistance_db/Mycobacterium_tuberculosis/mutations/mykrobe_annotated/correct.bed resistance_db/Mycobacterium_tuberculosis/mutations/miotto_high_moderate_minimum_confidence_annotated/correct.bed resistance_db/Mycobacterium_tuberculosis/mutations/walker_resistant_annotated/correct.bed

RUN cp ${pipeline_folder}/*.tsv . && cp ${pipeline_folder}/config.yaml .

RUN mkdir -p core_genomes/parsnp/Mycobacterium_tuberculosis/

RUN echo '' > links/Myco-10_S10_L001_R1.fastq.gz && echo '' > links/Myco-10_S10_L001_R2.fastq.gz && echo '' > core_genomes/parsnp/Mycobacterium_tuberculosis/parsnp.xmfa
#&& echo '' > ${main}/data/references/mash_sketch_human.msh

RUN wget -qO- https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh > ${main}/data/references/mash_sketch.msh

RUN mkdir -p ${main}/data/references/centrifuge_db/

RUN wget -qO- ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz | tar xvz -C ${main}/data/references/centrifuge_db/

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml report/multiqc_assembly/multiqc_report.html samples/M10/resistance/mykrobe.tsv report/figures/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_no_st.svg report/figures/gatk_gvcfs/core_parsnp_538048/bwa/distances_in_snp_mst_no_st.svg report/figures/gatk_gvcfs/full_genome_M10_assembled_genome/bwa/distances_in_snp_mst_no_st.svg virulence_summary.xlsx report/typing/mlst/summary.xlsx report/resistance/rgi_summary.xlsx report/resistance/mykrobe_summary.xlsx report/figures/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_no_st.svg report/figures/gatk_gvcfs/full_genome_538048/bwa/phylogeny_no_st.svg report/figures/gatk_gvcfs/core_parsnp_538048/bwa/phylogeny_no_st.svg report/contamination/mash/assembly/distances_formated.xlsx samples/M10/contamination/centrifuge/report.tsv -j 4 --config species="Mycobacterium_tuberculosis" && conda clean --all --yes

RUN echo '' > links/ERR2130394_fastq.gz
RUN echo '' > links/SRR6936587_fastq.gz

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml epidemiology --config species="Mycobacterium_tuberculosis" && conda clean --all --yes

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml resistance --config species="Mycobacterium_tuberculosis" && conda clean --all --yes

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml virulence --config species="Mycobacterium_tuberculosis" && conda clean --all --yes

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml epidemiology --config species="Mycobacterium_tuberculosis" && conda clean --all --yes

RUN patch /opt/conda/envs/9d0a6ae9/lib/python2.7/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch

#RUN mkdir -p ${main}/data/references/1493941/

#RUN  /bin/bash -c 'source activate /opt/conda/envs/db3680fa/ && efetch -db assembly -id 1493941 -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | sed "s/\(\/GCF_.*\)/\\1\\1_genomic.fna.gz/" | xargs -I % wget -qO- % | gzip -d > ${main}/data/references/1493941/genome_fasta.fasta'

#RUN  /bin/bash -c 'source activate /opt/conda/envs/c327f08f/ && mash sketch -o ${main}/data/references/mash_sketch_human.msh ${main}/data/references/1493941/genome_fasta.fasta'

#RUN rm -rf ${main}/data/references/1493941/

RUN rm -rf links/ &&  rm core_genomes/parsnp/Mycobacterium_tuberculosis/parsnp.xmfa && rm config.yaml && rm *.tsv

RUN mkdir ${main}/data/analysis/

WORKDIR ${main}/data/analysis/

RUN wget https://card.mcmaster.ca/download/0/broadstreet-v3.0.0.tar.gz --no-check-certificate && tar xf broadstreet-v3.0.0.tar.gz

RUN /bin/bash -c 'source activate /opt/conda/envs/5181c089/ && rgi load -i card.json'

RUN chmod -R 777 /opt/conda/envs/5181c089/lib/python3.6/site-packages/app/_data

# give write acess to rgi database
RUN chmod -R 777 /opt/conda/envs/5181c089/lib/python3.6/site-packages/app/_db/

RUN chmod -R 777 /opt/conda/envs/9d0a6ae9/lib/python2.7/site-packages/mykatlas/

RUN chown -R pipeline_user ${main}/

USER pipeline_user

WORKDIR ${main}/data/analysis/
