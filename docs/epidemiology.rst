
.. _epidemiology:
   
Epidemiology
============

Depends on the :ref:`assembly_quality` workflow (for determining the Sequence Types).

----------
Parameters
----------

* ``minimum_coverage_for_calling``: minimum of coverage for considering a genomic position when counting differences between samples. Any position (whether the genotype is identical to the reference, ie `GT=0` in the vcf, or different, ie `GT=1`) having a lower coverage will be masked
* ``minimum_alternate_fraction_for_calling``: minimum ratio of observations favouring a SNP over observations not favouring a SNP. Any position (`GT=0` or `GT=1`) not meeting this criteria will also be masked  
.. * ``ref_ids_for``

--------------------
Available Genotypers
--------------------

Genotyping can be performed with two different softwares:
  
.. toctree::
   :glob:

   genotypers/*


---------      
Filtering   
---------

.. toctree::
   :glob:
   
   filtering/*


------------
Deliverables
------------

* ``typing/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/distance_snp_mst_no_st.svg``: Minimum spanning tree of the distance in snps between every sample over the core genome as defined by ridom or enterobase. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``typing/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/distance_snp_mst_with_st.svg`` can be generated with the ST of each sample.
* ``phylogeny/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/phylogeny_no_st.svg``: A phylogeny based on the alignments of the core SNPs, using RAxML. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``phylogeny/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/phylogeny_with_st.svg`` can be generated with the ST of each sample.
  
* ``quality/multiqc/mapping_to_{reference_genome}/multiqc_report.html``: multiqc report of **qualimap**, **fastqc** and **trimmomatic** of every samples when mapping against the reference. Check for quality control.

