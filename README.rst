Routine procedures for diagnostic purposes using microbial genomics and metagenomics.

Workflows for epidemiology, anti-microbial resistance genotyping and virulence factors identification have been implemented using the `Snakemake <http://snakemake.readthedocs.io/en/stable/>`_ workflow management system with `bioconda <https://bioconda.github.io/>`_ integration for software dependency. `Docker <https://hub.docker.com/r/metagenlab/diag_pipelines/>`_ images of main releases are available.


Dependencies
============

------
Docker
------

Install the CE version following these `instructions <https://docs.docker.com/install/linux/docker-ce/ubuntu/>`_ for ubuntu. Also make sure you have created the docker group and that you can run docker without sudo following these `instruction <https://docs.docker.com/install/linux/linux-postinstall/>`_. If you can't have access to the internet when inside a Docker container, apply those `changes <https://docs.docker.com/install/linux/linux-postinstall/#disable-dnsmasq>`_.

.. code-block:: bash

   docker run hello-world
   docker pull metagenlab/diag_pipelines:latest
   docker run -t --rm metagenlab/diag_pipelines:latest sh -c "ping www.google.com"

Our Docker image is fit for a user called ``pipeline_user`` whose UID is ``1080``. It is advised to create this user on your computer before using the Docker image to run your analysis.


.. code-block:: bash

   sudo useradd -G docker,sudo -u 1080 pipeline_user
   sudo mkdir /home/pipeline_user/
   sudo chown pipeline_user -R /home/pipeline_user/
   sudo passwd pipeline_user

Alternatively, you can run the Docker as root (``--user root``) but the created folders will belong to the root user of your computer.

General use
===========
Once you have pulled the docker image on your computer:

.. code-block:: bash

   docker run -t --rm \
   --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind \
   metagenlab/diag_pipelines:latest \
   sh -c 'snakemake --snakefile $pipeline_folder/workflows/full_pipeline.rules \
   --use-conda --conda-prefix $conda_folder --configfile config.yaml'

Update the config file for your needs. If you have read files you want to analyse, they should be stored in the ``links`` folder from your current working directory.

Run one of the 4 workflows
============================

The pipeline implement four main workflows.

1. `Epidemiological analysis`_
2. `Annotation of virulence factors`_
3. `Annotation of resistance markers`_
4. `Characterization of one or multiple strains`_

An *html* report summarizing results is generated upon completion of the workflow.

------
Epidemiological analysis
------

.. code-block:: bash

		docker run -t --rm \
		--mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind \
		metagenlab/diag_pipelines:latest \
		sh -c 'snakemake --snakefile $pipeline_folder/workflows/full_pipeline.rules\
		--use-conda --conda-prefix $conda_folder --configfile config.yaml\
		epidemiology'

This will perform quality checks, map reads against one or multiple reference genome, calculate pairwise number of SNPs
and generate a minimum spanning tree. The reference genome can be one of the assembled genome or an assembly available on
the NCBI website. The analysis can also be restricted to the core genome as defined by existing cgMLST schemes or by
computing a custom core genome with help of parsnp_ (see documentation https://metagenlabdiag-pipelines.readthedocs.io/en/latest/core_genomes.html).

------
Annotation of virulence factors
------

.. code-block:: bash

		docker run -t --rm \
		--mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind \
		metagenlab/diag_pipelines:latest \
		sh -c 'snakemake --snakefile $pipeline_folder/workflows/full_pipeline.rules\
		--use-conda --conda-prefix $conda_folder --configfile config.yaml\
		virulence'

This will perform quality checks, assemble the genome and search for known virulence factors from the VFDB_ database.

------
Annotation of resistance markers
------

.. code-block:: bash

		docker run -t --rm \
		--mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind \
		metagenlab/diag_pipelines:latest \
		sh -c 'snakemake --snakefile $pipeline_folder/workflows/full_pipeline.rules\
		--use-conda --conda-prefix $conda_folder --configfile config.yaml\
		resistance'

This will perform quality checks, assemble the genome and search for known antibiotic resistance determinants with help of the `rgi software`_ and CARD_ database.

------
Characterization of one or multiple strains
------

.. code-block:: bash

		docker run -t --rm \
		--mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind \
		metagenlab/diag_pipelines:latest \
		sh -c 'snakemake --snakefile $pipeline_folder/workflows/full_pipeline.rules\
		--use-conda --conda-prefix $conda_folder --configfile config.yaml\
		strain_characterization'

This will perform quality checks, assemble the genome and search for known antibiotic resistance determinants with
help of the `rgi software`_ and CARD_ database and search for known virulence factors from the VFDB_ database.



Generating specific files of interest
============================

If you want to execute a specific analysis, you can request files of interest for a particular analysis. Consult the
full documentation to know what files can be generated (http://metagenlabdiag-pipelines.readthedocs.io/en/latest/ ).
Main examples are provided below:

.. code-block:: bash

   docker run -t --rm \
   --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind \
   metagenlab/diag_pipelines:latest \
   sh -c 'snakemake --snakefile $pipeline_folder/workflows/full_pipeline.rules\
   --use-conda --conda-prefix $conda_folder --configfile config.yaml\
   report/multiqc_assembly/multiqc_report.html'

This will assemble and annotate every samples, and generate a multiqc report for all samples.

.. code-block:: bash

   docker run -t --rm \
   --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind \
   metagenlab/diag_pipelines:latest \
   sh -c 'snakemake --snakefile $pipeline_folder/workflows/virulence.rules\
   --use-conda --conda-prefix $conda_folder --configfile config.yaml\
   virulence_summary.xlsx'

This will generate a summary excel file for the virulence factors of the samples, based on the virulence factors annotated in the file defined on the config file.

.. code-block:: bash

   docker run -t --rm \
   --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind \
   metagenlab/diag_pipelines:latest \
   sh -c 'snakemake --snakefile $pipeline_folder/workflows/typing.rules\
   --use-conda --conda-prefix $conda_folder --configfile config.yaml\
   typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp.xlsx'

This will generate a snp-distance matrix of all samples, only on the core genome defined by ridom of the species defined in the `species` variable of the config file, mapped with bwa on the reference genome used by ridom (which is *Staphylococcus aureus* COL substrain, `id` 33148 from the `NCBI Assembly database <https:/www.ncbi.nlm.nih.gov/assembly/>`_).

.. code-block:: bash

   docker run -t --rm \
   --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind \
   metagenlab/diag_pipelines:latest \
   sh -c 'snakemake --snakefile $pipeline_folder/workflows/resistance.rules\
   --use-conda --conda-prefix $conda_folder --configfile config.yaml\
   report/typing/mlst/summary.xlsx'

This will generate an Excel summary file of the MLST of all samples, based on the software `mlst <https://github.com/tseemann/mlst>`_.


----------------
All Deliverables
----------------

Here is a list of all deliverables currently available:


.. _VFDB:  http://www.mgc.ac.cn/VFs/
.. _rgi software: https://card.mcmaster.ca/analyze/rgi
.. _CARD: https://card.mcmaster.ca/
.. _parsnp: https://harvest.readthedocs.io/en/latest/content/parsnp.html
