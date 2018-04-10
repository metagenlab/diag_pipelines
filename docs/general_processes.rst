.. general_processes:
   
============
Generalities
============

------------------
Defining variables
------------------


As a general rules, any ``variable`` referenced in this documentation must be either:

* Defined in the yaml config file that is passed to snakemake by ``--configfile``
* Defined directly in the snakemake command by ``--config variable=$value``

General variables
-----------------

* Pass ``--cores number_of_cores`` to the snakemake command to define the number of cores you want to use for your analysis

An example of each variable value is provided in ``config.yaml`` of the github repository.

  
-----------------
Logging functions
-----------------

Logging processes are defined in the file :file:`workflows/logging.rules`.

Parameters
----------

The variable ``logging_folder`` can be defined in the ``config.yaml`` or passed to snakemake with ``--config``. If it is not defined, the default value will be the folder ``logs/``.

Saved files
-----------

Each time an effective snakemake run is started, a folder named with the current UTC datetime is created. A variable number of files will be copied there, so that replication of the run is possible:

* The snakefile passed to snakemake
* The config file
* The full command used, copied into the file ``cmd.txt``
* The parameter files defining the SRA and the local samples, if they exist
 
The logs of every command run during the execution of the workflow will then be stored in this folder.
  
------------------------
Determining sample names
------------------------

Sample naming and matching to fastq files are handled in the file :file:`workflows/making_sample_dataset.rules`.

Parameters
----------

* The variable ``link_directory`` can be defined. Its default value is ``links/``. Read data for local samples will be searched in this directory.
* :ref:`local_samples` can be defined based on a tabulated file whose full path can be passed to the variable ``local_samples``
* :ref:`sra_samples` can be defined based on the tabulated file whose full path can be passed to the variable ``sra_samples``
  

.. _local_samples:
     
Local samples
-------------

The tabulated file whose path is defined by ``local_samples`` must contain at least two columns: `SampleName` and `ScientificName`.

.. csv-table:: Local data example
   :header: "SampleName", "ScientificName"
   
   "S10","Staphylococcus aureus"
   "S1","Staphylococcus aureus"


For each entry, there must be in the folder defined by the ``link_directory`` variable, two files (for paired reads) or only one (for single reads) whose filename starts by one and only one entry of the `SampleName` columns. For instance, the files ``S10_001_R1_L001.fastq.gz`` and ``S10_001_R2_L001.fastq.gz`` in the folder defined by the ``link_directory`` variable will be matched to the sample name ``S10``. The matching is performed by using regular expressions to end the search at non alphanumeric characters or by the end of the word, thus the sample name ``S1`` will actually not match ``S10_001_R1_L001.fastq.gz`` nor ``S10_001_R2_L001.fastq.gz``.

If needed, an `OldSampleName` column can be added to the file, when the read filenames and the desired new sample names can not be matched simply by testing the identity at the start of both names. 

.. csv-table:: Local data example with old sample names
   :header: "SampleName", "ScientificName", "OldSampleName"
   
   "S10","Staphylococcus aureus","Staaur-10"
   "S1","Staphylococcus aureus","Staaur-1"	
      
In this case, the files ``Staaur-10_S10_L001_R1_001.fastq.gz`` and ``Staaur-10_S10_L001_R2_001.fastq.gz`` in the folder defined in ``link_directory`` will be matched to the sample name ``S10``. Similarly, ``Staaur-1`` will actually not match ``Staaur-10_S10_L001_R1_001.fastq.gz``.

.. _sra_samples:

SRA samples
-----------

The tabulated file whose path is defined by ``sra_samples`` can be the RunInfo files that are downloadable directly through the `SRA NCBI <https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA295367>`_ and can be passed without any modification. If the variable ``use_library_name`` with any value is passed during execution, the column `LibraryName` will be used for naming the samples, instead of `SampleName` (which can be useful for badly formatted Sra Run Info files). If you want to format your own SRA samples without a RunInfo file from the NCBI, 4 columns must be defined:

.. csv-table:: SRA data example
   :header: "Run","SampleName", "LibraryLayout", "ScientificName"

   "ERR2130394","SE1","single","Mycobacterium tuberculosis"
   "SRR6936587","XTB-14-012","single","Mycobacterium tuberculosis"


.. _workflows:

-------------------
Workflows and Rules
-------------------

Current available workflows are implemented in the folder ``workflows``. Each workflow will depend on ``rules``, stored in the folder of the same name, and can also depend on other workflows. ``rules`` are sorted with respect to their general function in different folders.

.. toctree::

   all_rules.rst

