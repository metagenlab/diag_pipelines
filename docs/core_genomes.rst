
.. _core_genome:

=========================
Core genome determination
=========================

Core genomes can be calculated by three different means.


-----
Ridom
-----

cgMLST scheme from `ridom <http://www.cgmlst.org/ncs>`_ can be extracted directly for theses species
  
.. csv-table:: Available cgMLST schemes from ridom
   :header: "Species", "Taxonomy ID", "Ridom ID", "Reference genome assembly ID"

   "*Staphylococcus_aureus*","1280","141106","33148"
   "*Mycobacterium_tuberculosis*","1773","741110","538048"
   "*Listeria_monocytogenes*","1639","690488","264498"
   "*Escherichia_coli*","562","5064703","79781"
   "*Klebsiella_pneumoniae*","573","2187931","31388"
   "*Enterococcus_faecium*","1352","991893","526908"
   "*Acinetobacter_baumannii*","470","3956907","39528"
   "*Legionella_pneumophila*","446","1025099","30068"

A bed file is constructed from the locus target file, using coordinates from the start and length columns of the csv file file available on the `ridom website <http://www.cgmlst.org/ncs/schema/3956907/locus/?content-type=csv>`_. 

Example
-------

.. code-block:: bash

   snakemake --snakefile $pipeline_folder/workflows/core_genome/make_ridom.rules --use-conda --conda-prefix $conda_folder --config species="Staphylococcus_aureus" -f all

will create a BED file in ``core_genomes/cgMLST/Staphylococcus_aureus.bed`` which defines the core genomic regions in the genome of the assembly ID 33148 (*Staphylococcus aureus* COL). 

----------
Enterobase
----------

cgMLST scheme from `enterobase <http://enterobase.warwick.ac.uk/>`_ is extracted for *Salmonella enterica*:



.. csv-table:: Available cgMLST schemes from enterobase
   :header: "Species", "Taxonomy ID", "Enterobase ID", "Reference genome assembly ID", "Scheme"

   "*Salmonella_enterica*","28901","SALwgMLST","359488","cgMLSTv1"

A bed file for the reference genome `359488 <https://www.ncbi.nlm.nih.gov/assembly/GCF_000027025.1/>`_, based on the locus tag present in this genome is constructed. For instance, over the 3002 loci of the *Salmonella* cgMLSTv1, 69 come from a different genome than the reference 359488.

Example
-------

.. code-block:: bash

   snakemake --snakefile $pipeline_folder/workflows/core_genome/make_enterobase.rules --use-conda --conda-prefix $conda_folder --config species="Salmonella enterica" -f all

will create a BED file in ``core_genomes/cgMLST/Salmonella_enterica.bed`` defining the core genomic regions in the genome of the assembly ID 359488 (*Salmonella enterica* subsp. enterica serovar Typhimurium str. D23580).
   

------   
ParSNP
------

For species unavailable on either resource, core genome can be calculated using parsnp and the complete genomes of the species available on RefSeq. As ParSNP is not available on bioconda, the binary must be downloaded from the `ParSNP website <http://harvest.readthedocs.io/en/latest/content/parsnp/quickstart.html>`_ and placed in your $PATH. 

Example
-------

.. code-block:: bash
		
   snakemake --config species="Morganella morganii" taxid="582" --snakefile $pipeline_folder/workflows/core_genomes/make_parsnp.rules --use-conda --conda-prefix $conda_folder -f all

will calculate the core genome with parSNP with every complete genome of *Morganella morganii* available in `RefSeq <https://www.ncbi.nlm.nih.gov/refseq/>`_. The ``taxid`` value must be the `taxonomy ID <https://www.ncbi.nlm.nih.gov/taxonomy/>`_ for the species defined. The resulting file will be located in ``core_genomes/parsnp/Morganella_morganii/parsnp.xmfa``.
