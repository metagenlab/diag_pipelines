.. _workflows:

=========
Workflows
=========

Current available workflows are implemented in the folder ``workflows``. Each workflow will depend on ``rules``, stored in the folder of the same name, and can also depend on other workflows. ``rules`` are sorted with respect to their general function in different folders.


.. _core_genome:

Core genome determination
=========================

Core genomes can be calculated by three different means.


-----
Ridom
-----

cgMLST scheme from `ridom <http://www.cgmlst.org/ncs>`_ can be extracted directly for theses species
  
.. csv-table:: Available cgMLST schemes from ridom
   :header: "Species", "Taxonomy ID", "Ridom ID", "Reference genome assembly ID"

   "*Staphylococcus aureus*","1280","141106","33148"
   "*Mycobacterium tuberculosis*","1773","741110","538048"
   "*Listeria monocytogenes*","1639","690488","264498"
   "*Escherichia coli*","562","5064703","79781"
   "*Klebsiella pneumoniae*","573","2187931","31388"
   "*Enterococcus faecium*","1352","991893","526908"
   "*Acinetobacter baumannii*","470","3956907","39528"
   "*Legionella pneumophila*","446","1025099","30068"

A bed file is constructed from the locus target file, constructing coordinates from the start and length columns of the csv file file available on the `ridom website <http://www.cgmlst.org/ncs/schema/3956907/locus/?content-type=csv>`_. 

Example
-------

.. code-block:: console

   snakemake --snakefile $pipeline_folder/workflows/core_genome/make_ridom.rules --use-conda --conda-prefix $conda_folder --config species="Staphylococcus aureus" -f all

will create a BED file in ``core_genomes/Staphylococcus_aureus/ridom/33148.bed`` which defines the core genomic regions in the genome of the assembly ID 33148 (*Staphylococcus aureus* COL). 

----------
Enterobase
----------

cgMLST scheme from `enterobase <http://enterobase.warwick.ac.uk/>`_ is extracted for *Salmonella enterica*:



.. csv-table:: Available cgMLST schemes from enterobase
   :header: "Species", "Taxonomy ID", "Enterobase ID", "Reference genome assembly ID", "Scheme"

   "*Salmonella enterica*","28901","SALwgMLST","359488","cgMLSTv1"

A bed file for the reference genome `359488 <https://www.ncbi.nlm.nih.gov/assembly/GCF_000027025.1/>`_, based on the locus tag present in this genome is constructed. For instance, over the 3002 locus of the *Salmonella* cgMLSTv1, 69 come from a different genome than the reference 359488.

Example
-------

.. code-block:: console

   snakemake --snakefile $pipeline_folder/workflows/core_genome/make_enterobase.rules --use-conda --conda-prefix $conda_folder --config species="Salmonella enterica" -f all

will create a BED file in ``core_genomes/Salmonella_enterica/enterobase/359488.bed`` defining the core genomic regions in the genome of the assembly ID 359488 (*Salmonella enterica* subsp. enterica serovar Typhimurium str. D23580).
   

------   
ParSNP
------

For species unavailable on either resource, core genome can be calculated using parsnp and the complete genomes of the species available on RefSeq. As ParSNP is not available on bioconda, the binary must be downloaded from the `ParSNP website <http://harvest.readthedocs.io/en/latest/content/parsnp/quickstart.html>`_ and placed in your $PATH. 

Example
-------

.. code-block:: bash
		
   snakemake --config species="Morganella morganii" taxid="582" --snakefile $pipeline_folder/workflows/core_genomes/make_parsnp.rules --use-conda --conda-prefix $conda_folder -f all

will calculate the core genome with parSNP with every complete genome of *Morganella morganii* available in `RefSeq <https://www.ncbi.nlm.nih.gov/refseq/>`_. The ``taxid`` value must be the `taxonomy ID <https://www.ncbi.nlm.nih.gov/taxonomy/>`_ for the species defined. The resulting file will be located in ``core_genomes/Morganella_morganii/parsnp/parsnp.xmfa``.

.. _assembly_quality:
     
Assembly and quality
====================


Aggregates rules for assembling genomes and performing various quality control checks.

----------
Parameters
----------

* ``cov_cutoff``: contigs whose coverage is below this cutoff will be excluded from the final assembly
* ``adapter_file_name``: look for the adaptor for this library preparation kit (possible `values <https://github.com/timflutre/trimmomatic/tree/master/adapters>`_)
* ``adapter_removal_param1``, ``adapter_removal_param2``, ``adapter_removal_param3``: parameters for adapter trimming (`reference <http://www.usadellab.org/cms/index.php?page=trimmomatic>`_)
* ``minimum_quality_base``: leading and trailing bases below this quality will be removed
* ``minimum_read_length``: reads shorter than this threshold after trimming will be discarded (be careful when using reads from SRA!)

------------
Deliverables
------------
 
* ``quality/multiqc/self_genome/multiqc_report.html``: quality control report based on the results of **fastqc**, **trimmomatic**, **qualimap**, **quast** and **prokka** for every sample
* ``samples/{sample_name}/annotation/``: folder containing all annotation files from **prokka**

.. _resistance:

Resistance
==========

Depends on the :ref:`assembly_quality` workflow.

----------
Parameters
----------

* ``resistance_prediction_softwares``: list of software for genetic resistance assessment. Possible values: ``mykrobe`` and ``rgi``.
  
.. * ``currated_resistance_genes``: file of trusted genes involved in resistance. An example is available in the folder ``data/mycobacterium/db/``

------------
Deliverables
------------

* ``samples/{sample_name}/annotation/resistance/rgi.tsv``: results files for RGI 
* ``samples/{sample_name}/annotation/resistance/mykrobe.tsv``: results file for mykrobe


.. _virulence:
  
Virulence
=========

Depends on the :ref:`assembly_quality` workflow.

----------
Parameters
----------

* ``virulence_factors``: file with list of uniprot accession of virulence factors. An example is available in the folder ``data/staph/db/``

------------
Deliverables
------------

*  ``virulence_summary.xlsx``: summary of virulence proteins found in every samples.


.. _epidemiology:
   
Epidemiology
============

Depends on the :ref:`assembly_quality` workflow (for determining the Sequence Types).

----------
Parameters
----------

* ``minimum_coverage_for_calling``: minimum of coverage for considering a genomic position when counting differences between samples. Any position (SNP or non-SNP when compared to the reference) having a lower coverage will be masked
* ``minimum_alternate_fraction_for_calling``: minimum ratio of observations favouring a SNP over observations not favouring a SNP. Any SNPs not meeting this criteria will also be masked
  
.. * ``ref_ids_for``

---------------------
Available SNP callers
---------------------

SNPs can be called with two different softwares:

* :ref:`freebayes`
* :ref:`gatk`
   
------------
Deliverables
------------

* ``typing/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/distance_snp_mst_no_st.svg``: Minimum spanning tree of the distance in snps between every sample over the core genome as defined by ridom or enterobase. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``typing/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/distance_snp_mst_with_st.svg`` can be generated with the ST of each sample.
* ``phylogeny/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/phylogeny_no_st.svg``: A phylogeny based on the alignments of the core SNPs, using RAxML. Available species and values for reference genomes are listed in the files in ``data/core_genome_dbs/``. If the species under consideration has a multiple locus sequence type available, ``phylogeny/{snp_caller}/core_{ridom or enterobase}/{reference_genome}/bwa/phylogeny_with_st.svg`` can be generated with the ST of each sample.
  
* ``quality/multiqc/mapping_to_{reference_genome}/multiqc_report.html``: multiqc report of **qualimap**, **fastqc** and **trimmomatic** of every samples when mapping against the reference. Check for quality control.
  
.. toctree::
