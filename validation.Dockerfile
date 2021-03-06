FROM metagenlab/diag_pipelines:2.1.2

# TESTING RESISTANT MTB ISOLATES WITH ALL 3 METHODS

WORKDIR ${pipeline_folder}

RUN git pull

RUN git fetch

RUN mkdir -p ${main}/data/validation/MTB-XTR/

WORKDIR ${main}/data/validation/MTB-XTR/

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

RUN ln -s ${main}/data/resistance_db/ resistance_db

RUN vdb-config --restore-defaults

# THESE ARE SRAS I SELECTED FROM TB PORTALS NIAID FROM SAMPLES THAT ARE XTR

RUN /bin/bash -c 'source activate /opt/conda/envs/db3680fa/ && esearch -db sra -query "SRR1158874[ID] OR SRR1158923[ID] OR SRR1158907[ID] OR SRR1158898[ID]" | efetch -db sra -format runinfo | sed "s/,/\t/g" > MTB-XTR.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 report/resistance/mykrobe_summary.xlsx report/resistance/rgi_summary.xlsx samples/XTB13-114/resistance/bwa/mykrobe_annotated/mutations.vcf samples/XTB13-114/resistance/bwa/rgi_annotated_full_2_0_0/mutations.vcf samples/XTB13-137/resistance/bwa/mykrobe_annotated/mutations.vcf samples/XTB13-137/resistance/bwa/rgi_annotated_full_2_0_0/mutations.vcf samples/XTB13-134/resistance/bwa/mykrobe_annotated/mutations.vcf samples/XTB13-134/resistance/bwa/rgi_annotated_full_2_0_0/mutations.vcf samples/XTB13-122/resistance/bwa/mykrobe_annotated/mutations.vcf samples/XTB13-122/resistance/bwa/rgi_annotated_full_2_0_0/mutations.vcf report/contamination/mash/assembly/distances_formated.xlsx --config sra_samples="MTB-XTR.tsv" species="Mycobacterium_tuberculosis"

# STAPHYLOCOCCUS AUREUS

RUN mkdir -p ${main}/data/validation/PRJEB7847/

WORKDIR ${main}/data/validation/PRJEB7847/

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

RUN /bin/bash -c 'source activate /opt/conda/envs/db3680fa/ && esearch -db sra -query "PRJEB7847[BIOPROJECT] AND \"Staphylococcus aureus\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 5 > PRJEB7847_Staphylococcus-aureus.tsv'

#TESTING MYKROBE ON 4 STAPH AUREUS

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 report/resistance/mykrobe_summary.xlsx --config sra_samples=PRJEB7847_Staphylococcus-aureus.tsv species="Staphylococcus_aureus"

#TESTING TYPING ON 4 STAPH AUREUS
# No rule to produce typing/gatk_gvcfs/cgMLST/bwa/distances_in_snp_mst_with_st.svg
RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 report/figures/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_with_st.svg report/figures/gatk_gvcfs/cgMLST/bwa/distances_in_snp_mst_with_st.svg phylogeny/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_with_st.svg report/resistance/mykrobe_summary.xlsx report/resistance/rgi_summary.xlsx typing/freebayes_joint_genotyping/cgMLST/bwa/vcfs/SAMEA3355592_-_SAMEA3355596_in_snp.vcf report/figures/gatk_gvcfs/full_genome_SAMEA3355592_assembled_genome/bwa/distances_in_snp_mst_with_st.svg --config sra_samples=PRJEB7847_Staphylococcus-aureus.tsv species="Staphylococcus_aureus"

#TESTING TYPING OF SERATIA MARCESCENS ISOLATES. ANNOTATED AS SINGLE BUT ACTUALLY PAIRED READS, BUT SINGLE READS RULES ARE TESTED IF WE LEAVE IT LIKE THIS

RUN mkdir -p ${main}/data/validation/PRJNA396838

WORKDIR ${main}/data/validation/PRJNA396838

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

RUN /bin/bash -c 'source activate /opt/conda/envs/db3680fa/ && esearch -db sra -query "PRJNA396838[BIOPROJECT] AND \"Serratia marcescens\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 5 > PRJNA396838_Serratia-marcescens.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 report/resistance/rgi_summary.xlsx report/figures/freebayes_joint_genotyping/full_genome_89161/bwa/distances_in_snp_mst_no_st.svg --config sra_samples="PRJNA396838_Serratia-marcescens.tsv" species="Serratia_marcescens"

RUN mkdir -p ${main}/data/validation/PRJEB12011

WORKDIR ${main}/data/validation/PRJEB12011

RUN ln -s ${main}/data/core_genomes/ core_genomes

RUN ln -s ${main}/data/references/ references

RUN /bin/bash -c 'source activate /opt/conda/envs/db3680fa/ && esearch -db sra -query "PRJEB12011[BIOPROJECT] AND \"Mycobacterium tuberculosis complex\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 2 > PRJEB12011_Mycobacterium-tuberculosis.tsv'

#TESTING MYKROBE BRADLEY 2015 PANEL, RGI, ON 2 M. TB ISOLATES

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 report/resistance/rgi_summary.xlsx report/resistance/mykrobe_summary.xlsx --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis" mykrobe_panel="bradley-2015"

RUN rm -rf samples/

RUN rm -rf resistance/

#TESTING WALKER 2015 PANEL WITH MYKROBE

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 report/resistance/mykrobe_summary.xlsx --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis" mykrobe_panel="walker-2015"

RUN rm -rf resistance/

#TESTING BOTH PANELS AT THE SAME TIME WITH MYKROBE

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --notemp --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 report/resistance/mykrobe_summary.xlsx --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"

#TYPING ON 4 M. TB ISOLATES

RUN /bin/bash -c 'source activate /opt/conda/envs/db3680fa/ && esearch -db sra -query "PRJEB12011[BIOPROJECT] AND \"Mycobacterium tuberculosis complex\"[ORGANISM]" | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 4 > PRJEB12011_Mycobacterium-tuberculosis.tsv'

RUN snakemake --snakefile ${pipeline_folder}/workflows/resistance.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 report/figures/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_no_st.svg report/figures/gatk_gvcfs/cgMLST/bwa/distances_in_snp_mst_no_st.svg report/figures/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_no_st.svg report/resistance/mykrobe_summary.xlsx typing/freebayes_joint_genotyping/cgMLST/bwa/vcfs/SAMEA3697004_-_SAMEA3697005_in_snp.vcf report/figures/gatk_gvcfs/full_genome_SAMEA3697006_assembled_genome/bwa/distances_in_snp_mst_no_st.svg --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"


RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 epidemiology --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 resistance --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 virulence --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"

RUN snakemake --snakefile ${pipeline_folder}/workflows/full_pipeline.rules --use-conda --conda-prefix $conda_folder --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -j 4 strain_characterization --config sra_samples=PRJEB12011_Mycobacterium-tuberculosis.tsv species="Mycobacterium_tuberculosis"
