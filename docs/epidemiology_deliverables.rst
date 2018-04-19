
The wildcard ``{snp_caller}`` in the following result files can have two values: ``freebayes_joint_genotyping`` or ``gatk_gvcfs``.

* ``typing/{snp_caller}/cgMLST/bwa/distances_in_snp_mst_no_st.svg``: minimum spanning tree of the distance in snps between every sample over the core genome as defined by ridom or enterobase. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``typing/{snp_caller}/cgMLST/bwa/distances_in_snp_mst_with_st.svg`` can be generated with the ST of each sample.
  
* ``phylogeny/{snp_caller}/cgMLST/bwa/phylogeny_no_st.svg``: a phylogeny based on the alignments of the core SNPs, using RAxML. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``phylogeny/{snp_caller}/cgMLST/bwa/phylogeny_with_st.svg`` can be generated with the ST of each sample.
  
* If you do not want or cannot use core genome schemes, ``typing/{snp_caller}/full_genome_33148/bwa/distances_in_snp_mst_no_st.svg`` will show the minimum spanning tree over the full genome of the assembly ID ``33148`` (*S. aureus COL* genome from NCBI).

* If you want to genotype with mapping over one of your own sequenced sample, ``typing/{snp_caller}/full_genome_S10_assembled_genome/bwa/distances_in_snp_mst_no_st.svg`` will show the minimum spanning tree when mapping onto the sample called ``S10``.
