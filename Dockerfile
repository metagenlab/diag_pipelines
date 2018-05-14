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

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_ridom.rules --use-conda --conda-prefix ${conda_folder} core_genomes/cgMLST/Staphylococcus_aureus.bed core_genomes/cgMLST/Mycobacterium_tuberculosis.bed core_genomes/cgMLST/Listeria_monocytogenes.bed core_genomes/cgMLST/Klebsiella_pneumoniae.bed core_genomes/cgMLST/Enterococcus_faecium.bed core_genomes/cgMLST/Acinetobacter_baumannii.bed core_genomes/cgMLST/Legionella_pneumophila.bed core_genomes/cgMLST/Clostridioides_difficile.bed -j 4 

RUN snakemake --snakefile ${pipeline_folder}/workflows/core_genomes/make_enterobase.rules --use-conda --conda-prefix ${conda_folder} core_genomes/cgMLST/Salmonella_enterica.bed core_genomes/cgMLST/Escherichia_coli.bed -j 4

RUN cp ${pipeline_folder}/*.tsv . 

RUN cp ${pipeline_folder}/config.yaml .

RUN echo '' > links/Myco-10_S10_L001_R1.fastq.gz

RUN echo '' > links/Myco-10_S10_L001_R2.fastq.gz

RUN mkdir -p core_genomes/parsnp/Mycobacterium_tuberculosis/

RUN echo '' > core_genomes/parsnp/Mycobacterium_tuberculosis/parsnp.xmfa

RUN wget -qO- https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh > ${main}/data/references/mash_sketch.msh

RUN echo '' > ${main}/data/references/mash_sketch_human.msh

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --create-envs-only --conda-prefix ${conda_folder} --configfile config.yaml quality/multiqc/assembly/multiqc_report.html samples/M10/resistance/mykrobe.tsv typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/core_parsnp_538048/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/full_genome_M10_assembled_genome/bwa/distances_in_snp_mst_no_st.svg virulence_summary.xlsx typing/mlst/summary.xlsx resistance/rgi_summary.xlsx resistance/mykrobe_summary.xlsx phylogeny/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_no_st.svg phylogeny/gatk_gvcfs/full_genome_538048/bwa/phylogeny_no_st.svg phylogeny/gatk_gvcfs/core_parsnp_538048/bwa/phylogeny_no_st.svg contamination/distances_formated.xlsx resistance/annotations/currated_db_isoniazid/correct.bed -j 4 --config species="Mycobacterium_tuberculosis"

RUN conda clean --all

RUN patch /opt/conda/envs/356da27d/lib/python3.6/site-packages/mykatlas/typing/typer/presence.py < ${pipeline_folder}/patches/mykrobe.patch

RUN mkdir -p ${main}/data/references/1493941/

RUN  /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && efetch -db assembly -id 1493941 -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | sed "s/\(\/GCF_.*\)/\\1\\1_genomic.fna.gz/" | xargs -I % wget -qO- % | gzip -d > ${main}/data/references/1493941/genome_fasta.fasta'

RUN  /bin/bash -c 'source activate /opt/conda/envs/c327f08f/ && mash sketch -o ${main}/data/references/mash_sketch_human.msh ${main}/data/references/1493941/genome_fasta.fasta'

RUN rm -rf ${main}/data/references/1493941/

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

# TESTING RESISTANT MTB ISOLATES WITH ALL 3 METHODS

RUN mkdir -p ${main}/data/validation/MTB-XTR/

WORKDIR ${main}/data/validation/MTB-XTR/

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

# THESE ARE SRAS I SELECTED FROM TB PORTALS NIAID FROM SAMPLES THAT ARE XTR

RUN /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && esearch -db sra -query "SRR1158874[ID] OR SRR1158923[ID] OR SRR1158907[ID] OR SRR1158898[ID]" | efetch -db sra -format runinfo | sed "s/,/\t/g" > MTB-XTR.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 2 resistance/mykrobe_summary.xlsx resistance/rgi_summary.xlsx samples/XTB13-114/resistance/mykrobe_annotated/mutations.vcf samples/XTB13-114/resistance/rgi_annotated_full_2_0_0/mutations.vcf samples/XTB13-137/resistance/mykrobe_annotated/mutations.vcf samples/XTB13-137/resistance/rgi_annotated_full_2_0_0/mutations.vcf samples/XTB13-134/resistance/mykrobe_annotated/mutations.vcf samples/XTB13-134/resistance/rgi_annotated_full_2_0_0/mutations.vcf samples/XTB13-122/resistance/mykrobe_annotated/mutations.vcf samples/XTB13-122/resistance/rgi_annotated_full_2_0_0/mutations.vcf contamination/distances_formated.xlsx --config sra_samples="MTB-XTR.tsv" species="Mycobacterium_tuberculosis" 

# STAPHYLOCOCCUS AUREUS 

RUN mkdir -p ${main}/data/validation/PRJEB7847/

WORKDIR ${main}/data/validation/PRJEB7847/

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

RUN /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && esearch -db sra -query "PRJEB7847[BIOPROJECT] AND \"Staphylococcus aureus\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 5 > PRJEB7847_Staphylococcus-aureus.tsv'

#TESTING MYKROBE ON 4 STAPH AUREUS

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 2 resistance/mykrobe_summary.xlsx --config sra_samples=PRJEB7847_Staphylococcus-aureus.tsv species="Staphylococcus_aureus"

#TESTING TYPING ON 4 STAPH AUREUS

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 2 typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_with_st.svg typing/gatk_gvcfs/cgMLST/bwa/distances_in_snp_mst_with_st.svg phylogeny/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_with_st.svg resistance/mykrobe_summary.xlsx resistance/rgi_summary.xlsx typing/freebayes_joint_genotyping/cgMLST/bwa/vcfs/14035468_patient_1_Run2_-_14035482_patient_5_in_snp.vcf typing/gatk_gvcfs/full_genome_14035468_patient_1_Run2_assembled_genome/bwa/distances_in_snp_mst_with_st.svg --config sra_samples=PRJEB7847_Staphylococcus-aureus.tsv species="Staphylococcus_aureus"

#TESTING TYPING OF SERATIA MARCESCENS ISOLATES. ANNOTATED AS SINGLE BUT ACTUALLY PAIRED READS, BUT SINGLE READS RULES ARE TESTED IF WE LEAVE IT LIKE THIS 

RUN mkdir -p ${main}/data/validation/PRJNA396838

WORKDIR ${main}/data/validation/PRJNA396838

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

RUN /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && esearch -db sra -query "PRJNA396838[BIOPROJECT] AND \"Serratia marcescens\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 5 > PRJNA396838_Serratia-marcescens.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 2 resistance/rgi_summary.xlsx typing/freebayes_joint_genotyping/full_genome_89161/bwa/distances_in_snp_mst_no_st.svg --config sra_samples="PRJNA396838_Serratia-marcescens.tsv" species="Serratia_marcescens"

RUN mkdir -p ${main}/data/validation/PRJEB12011

WORKDIR ${main}/data/validation/PRJEB12011

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

RUN /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && esearch -db sra -query "PRJEB12011[BIOPROJECT] AND \"Mycobacterium tuberculosis complex\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 2 > PRJEB12011_Mycobacterium-tuberculosis.tsv'

#TESTING MYKROBE BRADLEY 2015 PANEL, RGI, ON 2 M. TB ISOLATES

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 2 resistance/rgi_summary.xlsx resistance/mykrobe_summary.xlsx --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis" mykrobe_panel="bradley-2015"

RUN rm -rf samples/

RUN rm -rf resistance/

#TESTING WALKER 2015 PANEL WITH MYKROBE

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 2 resistance/mykrobe_summary.xlsx --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis" mykrobe_panel="walker-2015"

RUN rm -rf resistance/

#TESTING BOTH PANELS AT THE SAME TIME WITH MYKROBE

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --notemp --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 2 resistance/mykrobe_summary.xlsx --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"

#TYPING ON 4 M. TB ISOLATES

RUN /bin/bash -c 'source activate /opt/conda/envs/618592fe/ && esearch -db sra -query "PRJEB12011[BIOPROJECT] AND \"Mycobacterium tuberculosis complex\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 5 > PRJEB12011_Mycobacterium-tuberculosis.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 2 typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/cgMLST/bwa/distances_in_snp_mst_no_st.svg phylogeny/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_no_st.svg resistance/mykrobe_summary.xlsx typing/freebayes_joint_genotyping/cgMLST/bwa/vcfs/SAMEA3697004_-_SAMEA3697005_in_snp.vcf typing/gatk_gvcfs/full_genome_SAMEA3697006_assembled_genome/bwa/distances_in_snp_mst_no_st.svg --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"

RUN conda clean --all

RUN chown -R pipeline_user ${main}/

USER pipeline_user

WORKDIR ${main}/data/analysis/


