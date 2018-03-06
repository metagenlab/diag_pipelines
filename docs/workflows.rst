.. workflows

---------
Workflows
---------

Current available workflows are implemented in the folder ``workflows``. Each workflow will depend on ``rules``, stored in the folder of the same name, and can also depend on other workflows. ``rules`` are sorted with respect to their general function in different folders.

``workflows`` for generating **core genomes** of species are also included. They can have three different origins:

* The cgMLST scheme of `ridom <http://www.cgmlst.org/ncs>`_
* The cgMLST scheme of `enterobase <http://enterobase.warwick.ac.uk/>`_
* For species unavailable on either resource, core genome can be calculated using parsnp and the complete genomes of the species available on RefSeq


.. _assembly_quality:
     
Assembly and quality
--------------------
Aggregates rules for assembling genomes and performing various quality control checks.
Required parameters:

* ``cov_cutoff``: contigs whose coverage is below this cutoff will be excluded from the final assembly
* ``adapter_file_name``: look for the adaptor for this library preparation kit (possible `values <https://github.com/timflutre/trimmomatic/tree/master/adapters>`_)
* ``adapter_removal_param1``, ``adapter_removal_param2``, ``adapter_removal_param3``: parameters for adapter trimming (`reference <http://www.usadellab.org/cms/index.php?page=trimmomatic>`_)
* ``minimum_quality_base``: leading and trailing bases below this quality will be removed
* ``minimum_read_length``: reads shorter than this threshold after trimming will be discarded (be careful when using reads from SRA!)

Deliverables:
 
* ``quality/multiqc/self_genome/multiqc_report.html``: quality control report based on the results of **fastqc**, **trimmomatic**, **qualimap**, **quast** and **prokka** for every sample
* ``samples/{sample_name}/annotation/``: folder containing all annotation files from the ``prokka`` software

.. _resistance:

Resistance
----------
Depends on the :ref:`assembly_quality` workflow.

Required parameters:

* ``resistance_prediction_softwares``: list of software for genetic resistance assessment. Possible values: ``mykrobe`` and ``rgi``.
  
.. * ``currated_resistance_genes``: file of trusted genes involved in resistance. An example is available in the folder ``data/mycobacterium/db/``
    
Deliverables:

* ``samples/{sample_name}/annotation/resistance/rgi.tsv``: results files for RGI 
* ``samples/{sample_name}/annotation/resistance/mykrobe.tsv``: results file for mykrobe


.. _virulence:
  
Virulence
---------
Depends on the :ref:`assembly_quality` workflow.

Required parameters:

* ``virulence_factors``: file with list of uniprot accession of virulence factors. An example is available in the folder ``data/staph/db/``
  
Deliverables:

*  ``virulence_summary.xlsx``: summary of virulence proteins found in every samples.


.. _epidemiology:
   
Epidemiology
------------
Depends on the :ref:`assembly_quality` workflow (for ST assessment).

Required parameters:

* ``minimum_coverage_for_calling``: minimum of coverage for considering a genomic position when counting differences between samples. Any position (SNP or non-SNP when compared to the reference) having a lower coverage will be masked
* ``minimum_alternate_fraction_for_calling``: minimum ratio of observations favouring a SNP over observations not favouring a SNP. Any SNPs not meeting this criteria will also be masked
  
.. * ``ref_ids_for``

Deliverables:

* ``typing/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/distance_snp_mst_no_st.svg``: Minimum spanning tree of the distance in snps between every sample over the core genome as defined by ridom or enterobase. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``typing/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/distance_snp_mst_with_st.svg`` can be generated with the ST of each sample.
* ``phylogeny/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/phylogeny_no_st.svg``: A phylogeny based on the alignments of the core SNPs, using RAxML. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``phylogeny/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/phylogeny_with_st.svg`` can be generated with the ST of each sample.
  
* ``quality/multiqc/mapping_to_{reference_genome}/multiqc_report.html``: multiqc report of **qualimap**, **fastqc** and **trimmomatic** of every samples when mapping against the reference. Check for quality control.
  
.. toctree::
