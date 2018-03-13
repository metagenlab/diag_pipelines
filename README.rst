Routine procedures for diagnostic purposes using microbial genomics and metagenomics.

Workflows for epidemiology, anti-microbial resistance genotyping and virulence factors identification have been implemented using the `Snakemake <http://snakemake.readthedocs.io/en/stable/>`_ workflow management system with `bioconda <https://bioconda.github.io/>`_ integration for software dependency. `Docker <https://hub.docker.com/r/metagenlab/diag_pipelines/>`_ images of main releases are available.


Dependencies
============
Docker  

.. code-block:: none
		
   docker pull metagenlab/diag_pipelines:ring_trial_v0.1.3


General use
===========
Once you have pulled the docker image on your computer: 

.. code-block:: none
		
  docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml

Update the config file for your needs.

Generating files of interest
============================

The pipeline works by asking the generation of the files of interest for a particular analysis. Consult the full documentation to know what files can be generated. Main examples are provided below: 

.. code-block:: none
		
  docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml config.yaml quality/multiqc/self_genome/multiqc_report.html'``

This will assemble and annotate every samples, and generate a multiqc report for all samples. 

.. code-block:: none
		
   docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml virulence_summary.xlsx'``

This will generate a summary excel file for the virulence factors of the samples, based on the virulence factors annotated in the file defined on the config file.

.. code-block:: none
		
   docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml typing/freebayes_joint_genotyping/core_ridom/33148/bwa/distances_snp.xlsx'``

This will generate a snp-distance matrix of all samples, only on the core genome calculated with parsnp and with all complete genomes of the species defined in the `taxid` variable of the config file, mapped on the assembly (from the `NCBI Assembly database <https:/www.ncbi.nlm.nih.gov/assembly/>`_) whose `id` is 33148 (*Staphylococcus aureus* COL substrain, reference genome for the ridom cgMLST scheme) with bwa. 

.. code-block:: none
		
   docker run -t --rm --mount source="$(pwd)",target=/home/pipeline_user/data,type=bind metagenlab/diag_pipelines:ring_trial_v0.1.3 sh -c 'snakemake --snakefile $pipeline_folder/workflows/ring_trial/pipeline.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml typing/mlst/summary.xlsx'``

This will generate an Excel summary file of the MLST of all samples, based on the software `mlst <https:/github.com/tseemann/mlst>`_)

