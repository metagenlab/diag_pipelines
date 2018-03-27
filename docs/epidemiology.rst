.. _epidemiology:
   
Epidemiology
============

Depends on the :ref:`assembly_quality` workflow (for determining the Sequence Types).
The genotyping results depend on the quality assessment performed on the mapping to the reference genomes, thus each time genotyping is performed, a Multiqc report is available in ``quality/multiqc/mapping_to_{ref}/multiqc_report.html`` and the contamination results in ``samples/{sample}/contamination/mash/distances_formated.xlsx`` for each sample.

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

* ``typing/{snp_caller}/cgMLST/bwa/distance_snp_mst_no_st.svg``: Minimum spanning tree of the distance in snps between every sample over the core genome as defined by ridom or enterobase. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``typing/{snp_caller}/cgMLST/bwa/distance_snp_mst_with_st.svg`` can be generated with the ST of each sample.
  
* ``phylogeny/{snp_caller}/cgMLST/bwa/phylogeny_no_st.svg``: A phylogeny based on the alignments of the core SNPs, using RAxML. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``phylogeny/{snp_caller}/cgMLST/bwa/phylogeny_with_st.svg`` can be generated with the ST of each sample.
  
* If you don't want or can't use core genome schemes, ``typing/{snp_caller}/full_genome/33148/bwa/distance_snp_mst_no_st.svg`` will show the minimum spanning tree over the full genome of the assembly ID ``33148`` (*S. aureus COL* genome from NCBI).

* If you want to genotyping with mapping over one of your own sequenced sample, ``typing/{snp_caller}/full_genome/S10_assembled_genome/bwa/distance_snp_mst_no_st.svg`` will show the minimum spanning tree when mapping onto the sample called ``S10``
