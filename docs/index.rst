.. diag_pipelines documentation master file, created by
   sphinx-quickstart on Wed Feb 28 12:10:42 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========================================================
Documentation for the genomics and metagenomics workflows
=========================================================

Routine procedures for diagnostic purposes using microbial genomics and metagenomics.

Workflows for epidemiology, anti-microbial resistance genotyping and virulence factors identification have been implemented using the `Snakemake <http://snakemake.readthedocs.io/en/stable/>`_ workflow management system with `bioconda <https://bioconda.github.io/>`_ integration for software dependency. `Docker <https://hub.docker.com/r/metagenlab/diag_pipelines/>`_ images of main releases are available.


---------
Workflows
---------

Current available workflows are implemented in the folder **workflows**. Each workflow will depend on **rules** and can also depend on other workflows. **Rules** are sorted with respect to their general function in different folders.

Higher level functions are also included in **workflows**:

* ``logging.rules`` archives all logs, commands, and configuration files used every time snakemake is run (excluding dry runs
* ``making_sample_dataset.rules`` determines which local samples and sra are going to be used on the run, based on the configuration files and the fastq read files present
   
**Workflows** for generating **core genomes** of species are also included. They can have three different origins:

* The cgMLST scheme of `ridom <http://www.cgmlst.org/ncs>`_
* The cgMLST scheme of `enterobase <http://enterobase.warwick.ac.uk/>`_
* For species unavailable on either resource, core genome can be calculated using parsnp and the complete genomes of the species available on RefSeq


Assembly and quality
--------------------
Aggregates rules for assembling genomes and performing various quality control checks.
Required parameters:

* ``cov_cutoff``: contigs whose coverage is below this cutoff will be excluded from the final assembly
* ``adapter_file_name``: look for the adaptor for this library preparation kit (possible `values <https://github.com/timflutre/trimmomatic/tree/master/adapters>`_)
* ``adapter_removal_param1``, ``adapter_removal_param2``, ``adapter_removal_param3``: parameters for adapter trimming (`reference <http://www.usadellab.org/cms/index.php?page=trimmomatic>`_)
* ``minimum_quality_base``: leading and trailing bases below this quality will be removed
* ``minimum_read_length``: reads shorter than this threshold after trimming will be discarded (be careful when using reads from SRA!)

Final deliverable:
 
``quality/multiqc/self_genome/multiqc_report.html``: quality control report based on the results of **fastqc**, **trimmomatic**, **qualimap**, **quast** and **prokka** for every sample

Resistance
----------
Depends on the **assembly and quality** workflow.
Required parameters:

* ``resistance_prediction_softwares``: list of software for genetic resistance assessment. Possible values: ``mykrobe`` and ``rgi``.
  
.. * ``currated_resistance_genes``: file of trusted genes involved in resistance. An example is available in the folder ``data/mycobacterium/db/``
    
Final deliverable: ``samples/{sample_name}/annotation/resistance/rgi.tsv`` and ``samples/{sample_name}/annotation/resistance/mykrobe.tsv``: results file for each sample from each prediction software


Virulence
---------
Depends on the **assembly and quality** workflow.
Required parameters:

* ``virulence_factors``: file with list of uniprot accession of virulence factors. An example is available in the folder ``data/staph/db/``
  
Final deliverable: ``virulence_summary.xlsx``: summary of virulence proteins found in every samples.


Epidemiology
------------
Depends on the **assembly and quality** workflow (for ST assessment).
Required parameters:

* ``species_taxid``: ID from the NCBI Taxonomy database of the species
  
.. * ``ref_ids_for``

Findal deliverable: ``typing/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/distance_snp_mst.pdf``: Minimum spanning tree of the distance in snps between every sample over the core genome as defined by ridom or enterobase. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``.
   
.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
