.. _core_genome:

=========================
Core genome determination
=========================

Core genomes can be calculated by three different means. Core genomes definition from :ref:`ridom` and :ref:`enterobase` are already included in the Docker image, in the folder ``/home/pipeline_user/core_genomes/cgMLST``, and the references they use in ``/home/pipeline_user/references/``.


.. _ridom:
-----
Ridom
-----

cgMLST scheme from `ridom <http://www.cgmlst.org/ncs>`_ can be extracted directly for theses species
  
.. csv-table:: Available cgMLST schemes from ridom
   :header: "Species", "Taxonomy ID", "Ridom ID", "Reference genome assembly ID"

   "*Acinetobacter_baumannii*","470","3956907","39528"
   "*Enterococcus_faecium*","1352","991893","526908"
   "*Klebsiella_pneumoniae*","573","2187931","31388"
   "*Listeria_monocytogenes*","1639","690488","264498"
   "*Legionella_pneumophila*","446","1025099","30068"
   "*Mycobacterium_tuberculosis*","1773","741110","538048"
   "*Staphylococcus_aureus*","1280","141106","33148"

A bed file is constructed from the locus target file, using coordinates from the start and length columns of the csv file file available on the `ridom website <http://www.cgmlst.org/ncs/schema/3956907/locus/?content-type=csv>`_. 

Example
-------

.. code-block:: bash

   snakemake --snakefile $pipeline_folder/workflows/core_genome/make_ridom.rules \
   core_genomes/cgMLST/Staphylococcus_aureus.bed \
   --use-conda --conda-prefix $conda_folder

will create a BED file in ``core_genomes/cgMLST/Staphylococcus_aureus.bed`` which defines the core genomic regions in the genome of the assembly ID ``33148`` (*Staphylococcus aureus* COL). 

.. _enterobase:

----------
Enterobase
----------

cgMLST scheme from `enterobase <http://enterobase.warwick.ac.uk/>`_ is available for:



.. csv-table:: Available cgMLST schemes from enterobase
   :header: "Species", "Taxonomy ID", "Enterobase ID", "Reference genome assembly ID", "Scheme"

   "*Escherichia_coli*","562","ESCwgMLST","79781","cgMLSTv1"
   "*Salmonella_enterica*","28901","SALwgMLST","359488","cgMLSTv1"


A bed file for each reference genome, based on the locus tags present in this genome, is constructed. For instance, over the 3002 loci of the *Salmonella* cgMLSTv1, 69 come from a different genome than the reference ``359488``. For *E. coli*, only 15 loci are missing for the reference assembly (``79781``), out of 2498.

Example
-------

.. code-block:: bash

   snakemake --snakefile $pipeline_folder/workflows/core_genome/make_enterobase.rules \
   core_genomes/cgMLST/Salmonella_enterica.bed \
   --use-conda --conda-prefix $conda_folder

will create a BED file in ``core_genomes/cgMLST/Salmonella_enterica.bed`` defining the core genomic regions in the genome of the assembly ID ``359488`` (*Salmonella enterica* subsp. enterica serovar Typhimurium str. D23580).
   

------   
ParSNP
------

For species unavailable on either resource, core genome can be calculated using parsnp and the complete genomes of the species available on RefSeq. As ParSNP is not available on bioconda, the binary must be downloaded from the `ParSNP website <http://harvest.readthedocs.io/en/latest/content/parsnp/quickstart.html>`_ and placed in your $PATH. 

Example
-------

.. code-block:: bash
		
   snakemake --snakefile $pipeline_folder/workflows/core_genomes/make_parsnp.rules \
   core_genome/parsnp/Morganella_morganii/parsnp.xmfa \
   --use-conda --conda-prefix $conda_folder 

will calculate the core genome with parSNP with every complete genome of *Morganella morganii* available in `RefSeq <https://www.ncbi.nlm.nih.gov/refseq/>`_.


If you wish to create a new parSNP core genome definition with the Docker image (that include the ``parsnp`` binary), do not link any ``references`` or ``core_genomes`` from your working directory.
