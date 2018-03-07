

Dependencies
============
Docker  

``docker pull metagenlab/diag_pipelines:ring_trial_v0.1.3``


General use
===========
Once you have pulled the docker image on your computer: 

``docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml'``

Update the config file for your needs.

Generating files of interest
============================

The pipeline works by asking the generation of the files of interest for a particular analysis. Consult the full documentation to know what files can be generated. Main examples are provided below:


``docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml config.yaml quality/multiqc/self_genome/multiqc_report.html'``

This will assemble and annotate every samples present in the `links` folder, and generate a multiqc report.


``docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml virulence_summary.xlsx'``

This will generate a summary excel file for the virulence factors of the strains in the `links` folder, based on the virulence factors annotated in the file defined on the config file.

``docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml typing/freebayes_joint_genotyping/core_ridom/33128/bwa/distances_from_merged_pairs_of_vcf.xlsx'``

This will generate a snp-distance matrix of all samples, only on the core genome calculated with parsnp and with all complete genomes of the species defined in the `taxid` variable of the config file, mapped on the assembly (from https:/www.ncbi.nlm.nih.gov/assembly/) whose `id` is 33128 (*Staphylococcus aureus* COL substrain, reference genome for the ridom cgMLST scheme) with bwa


``docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml typing/mlst/summary.xlsx'``

This will generate an Excel summary file of the MLST of all samples, based on the software mlst (https:/github.com/tseemann/mlst)

