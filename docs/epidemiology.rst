.. _epidemiology:
   
Epidemiology
============

Depends on the :ref:`assembly_quality` workflow (for determining the Sequence Types).
The genotyping results depend on the quality assessment performed on the mapping to the reference genomes, thus each time genotyping is performed, a Multiqc report is available in ``quality/multiqc/mapping_to_{ref}/multiqc_report.html`` and the contamination results in ``contamination/distances_formated.xlsx`` for each sample.

----------
Parameters
----------

* ``minimum_coverage_for_calling``: minimum of coverage for considering a genomic position when counting differences between samples. Any position (whether the genotype is identical to the reference, ie `GT=0` in the vcf, or different, ie `GT=1`) having a lower coverage will be masked.
* ``minimum_alternate_fraction_for_calling``: minimum ratio of observations favouring an alternative allele over observations favouring the reference allele. Any position (`GT=0` or `GT=1`) not meeting this criteria will also be masked.
* ``snp_threshold``: pairs of samples having less than this number of SNP differences will be linked in the final minimum spanning tree
* ``minimum_spanning_tree_size``: size of the minimum spanning tree image, default is ``10``
* ``phylogeny_image_size``: size of the phylogeny image, default is ``800``
* ``species``
  
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


-----------------------
Calculating differences
-----------------------
Once the filtering has been performed, differences in snps are calculated between samples and against the reference. 


------------
Deliverables
------------

The wildcard ``{snp_caller}`` in the following result files can have two values: ``freebayes_joint_genotyping`` or ``gatk_gvcfs``.

* ``typing/{snp_caller}/cgMLST/bwa/distances_snp_mst_no_st.svg``: minimum spanning tree of the distance in snps between every sample over the core genome as defined by ridom or enterobase. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``typing/{snp_caller}/cgMLST/bwa/distances_snp_mst_with_st.svg`` can be generated with the ST of each sample.
  
* ``phylogeny/{snp_caller}/cgMLST/bwa/phylogeny_no_st.svg``: a phylogeny based on the alignments of the core SNPs, using RAxML. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``phylogeny/{snp_caller}/cgMLST/bwa/phylogeny_with_st.svg`` can be generated with the ST of each sample.
  
* If you do not want or cannot use core genome schemes, ``typing/{snp_caller}/full_genome/33148/bwa/distances_snp_mst_no_st.svg`` will show the minimum spanning tree over the full genome of the assembly ID ``33148`` (*S. aureus COL* genome from NCBI).

* If you want to genotype with mapping over one of your own sequenced sample, ``typing/{snp_caller}/full_genome/S10_assembled_genome/bwa/distances_snp_mst_no_st.svg`` will show the minimum spanning tree when mapping onto the sample called ``S10``.
