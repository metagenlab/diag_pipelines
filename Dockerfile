FROM continuumio/miniconda3:4.3.27

RUN conda config --add channels defaults

RUN conda config --add channels conda-forge

RUN conda config --add channels bioconda

RUN useradd -r -u 1080 pipeline_user

RUN apt-get install -y fontconfig

RUN conda install snakemake=4.8.0=py36_0 unzip

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

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} core_genomes/cgMLST/Staphylococcus_aureus.bed core_genomes/cgMLST/Mycobacterium_tuberculosis.bed core_genomes/cgMLST/Listeria_monocytogenes.bed core_genomes/cgMLST/Klebsiella_pneumoniae.bed core_genomes/cgMLST/Enterococcus_faecium.bed core_genomes/cgMLST/Acinetobacter_baumannii.bed core_genomes/cgMLST/Legionella_pneumophila.bed -j 4 

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_enterobase.rules --use-conda --conda-prefix ${conda_folder} core_genomes/cgMLST/Salmonella_enterica.bed core_genomes/cgMLST/Escherichia_coli.bed -j 4

RUN cp ${pipeline_folder}/*.tsv . 

RUN cp ${pipeline_folder}/config.yaml .

RUN echo '' > links/Myco-10_S10_L001_R1.fastq.gz

RUN echo '' > links/Myco-10_S10_L001_R2.fastq.gz

RUN mkdir -p core_genomes/parsnp/Mycobacterium_tuberculosis/

RUN echo '' > core_genomes/parsnp/Mycobacterium_tuberculosis/parsnp.xmfa

RUN wget -qO- https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh > ${main}/data/references/mash_sketch.msh

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml  quality/multiqc/assembly/multiqc_report.html samples/M10/resistance/mykrobe.tsv typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/core_parsnp_538048/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/full_genome_M10_assembled_genome/bwa/distances_in_snp_mst_no_st.svg virulence_summary.xlsx typing/mlst/summary.xlsx samples/M10/resistance/m_tuberculosis_resistance_genes_4_db_mutations/summary.xlsx resistance/rgi_summary.xlsx resistance/mykrobe_summary.xlsx phylogeny/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_no_st.svg phylogeny/gatk_gvcfs/full_genome_538048/bwa/phylogeny_no_st.svg phylogeny/gatk_gvcfs/core_parsnp_538048/bwa/phylogeny_no_st.svg -j 4

RUN patch /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch

RUN rm -rf links/

RUN rm core_genomes/parsnp/Mycobacterium_tuberculosis/parsnp.xmfa

RUN rm config.yaml

RUN rm *.tsv

RUN mkdir ${main}/data/analysis/

WORKDIR ${main}/data/analysis/

RUN wget https://card.mcmaster.ca/download/0/broadstreet-v1.1.9.tar.bz2 --no-check-certificate

RUN tar xf broadstreet-v1.1.9.tar.bz2

RUN /bin/bash -c 'source activate /opt/conda/envs/803617ab/ && rgi_load -i card.json'

RUN chmod -R 777 /opt/conda/envs/803617ab/lib/python2.7/site-packages/package/rgi/_db/

RUN chmod -R 777 /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/

RUN chown -R pipeline_user ${main}/

USER pipeline_user

RUN  /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && efetch -db assembly -id 1493941 -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | sed "s/\(\/GCF_.*\)/\\1\\1_genomic.fna.gz/" | xargs -I % wget -qO- % | gzip -d > ${main}/references/1493941/genome_fasta.fasta' 

RUN mkdir -p ${main}/data/validation/PRJEB7847/

WORKDIR ${main}/data/validation/PRJEB7847/

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

#NOW RUNNING TYPING, MLST, RESISTANCE (RGI+MYKROBE) FOR S. AUREUS

RUN /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && esearch -db sra -query "PRJEB7847[BIOPROJECT] AND \"Staphylococcus aureus\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 2 > PRJEB7847_Staphylococcus-aureus.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 resistance/mykrobe_summary.xlsx resistance/rgi_summary.xlsx virulence_summary.xlsx --config sra_samples=PRJEB7847_Staphylococcus-aureus.tsv species="Staphylococcus_aureus" virulence_factors=${pipeline_folder}/data/staph/db/virulence_factors.tsv virulence_percentage_identity_cutoff=80 virulence_coverage_cutoff=70

RUN /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && esearch -db sra -query "PRJEB7847[BIOPROJECT] AND \"Staphylococcus aureus\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 5 > PRJEB7847_Staphylococcus-aureus.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/typing.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_with_st.svg typing/gatk_gvcfs/cgMLST/bwa/distances_in_snp_mst_with_st.svg phylogeny/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_with_st.svg typing/freebayes_joint_genotyping/cgMLST/bwa/vcfs/14035468_patient_1_Run2_-_14035482_patient_5_in_snp.vcf typing/gatk_gvcfs/full_genome_14035468_patient_1_Run2_assembled_genome/bwa/distances_in_snp_mst_with_st.svg --config sra_samples=PRJEB7847_Staphylococcus-aureus.tsv species="Staphylococcus_aureus"

#NOW RUNNING TYPING, MLST, RESISTANCE (RGI+MYKROBE+HOMEMADE) FOR M. TUBERCULOSIS

RUN mkdir -p ${main}/data/validation/PRJEB12011

WORKDIR ${main}/data/validation/PRJEB12011

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

RUN /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && esearch -db sra -query "PRJEB12011[BIOPROJECT] AND \"Mycobacterium tuberculosis complex\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 2 > PRJEB12011_Mycobacterium-tuberculosis.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 resistance/rgi_summary.xlsx resistance/mykrobe_summary.xlsx resistance/m_tuberculosis_resistance_genes_4_db_mutations_summary.xlsx --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"

RUN rm -rf samples/

RUN rm -rf resistance/

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 resistance/mykrobe_summary.xlsx --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis" mykrobe_panel="walker-2015"

RUN /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && esearch -db sra -query "PRJEB12011[BIOPROJECT] AND \"Mycobacterium tuberculosis complex\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 5 > PRJEB12011_Mycobacterium-tuberculosis.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/typing.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/cgMLST/bwa/distances_in_snp_mst_no_st.svg phylogeny/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_no_st.svg typing/freebayes_joint_genotyping/cgMLST/bwa/vcfs/19-M-2-RUM_-_1633-10_in_snp.vcf typing/gatk_gvcfs/full_genome_1633-10_assembled_genome/bwa/distances_in_snp_mst_no_st.svg --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"

#NOW RUNNING SINGLE READS RULES

RUN mkdir -p ${main}/data/validation/single_reads/links

WORKDIR ${main}/data/validation/single_reads

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 resistance/mykrobe_summary.xlsx quality/multiqc/assembly/multiqc_report.html contamination/distances_formated.xlsx --config sra_samples=${pipeline_folder}/example_sra_samples.tsv species="Mycobacterium_tuberculosis"

WORKDIR ${main}/data/analysis/

RUN rm -rf *

RUN rm -rf ${main}/data/validation/


