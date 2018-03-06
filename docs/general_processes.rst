.. general_processes

-----------------
Logging functions
-----------------

Archiving processes are defined in the file :file:`../workflows/logging.rules`. The variable ``logging_folder`` must be defined in the ``config.yaml`` or passed to snakemake with ``--config``. Each time an effective snakemake run is started, a folder named with the current UTC datetime. A different number of files will be copied there, so that replication of the run is possible:

* The snakefile passed to snakemake
* The config file
* The full command used, copied into the file ``cmd.txt``
* The parameter files defined the SRA and the local samples, if they exist
 

------------------------
Determining sample names
------------------------

Samples for the run will be determined in the file :file:`../workflows/making_sample_dataset.rules`.


Local samples
-------------

Local samples will be determined based on the tabulated file whose full path must be passed to the variable ``local_samples`` in the ``config.yaml`` or through ``--config`` on the snakemake command. It must contain at least two columns: `SampleName` and `ScientificName`.

.. csv-table:: Local data example
   :header: "SampleName", "ScientificName"
   
   "S10","Staphylococcus aureus"
   "S1","Staphylococcus aureus"
   "S2","Staphylococcus aureus"
   "S3","Staphylococcus aureus"
   "S4","Staphylococcus aureus"
   "S5","Staphylococcus aureus"
   "S6","Staphylococcus aureus"
   "S7","Staphylococcus aureus"
   "S8","Staphylococcus aureus"
   "S9","Staphylococcus aureus"

For each entry, there must be in the folder defined in variable ``link_directory``, two files (for paired reads) or a single one (for single reads) whose filename starts with one and only one entry in the `SampleName` columns. For instance, the files ``S10_001_R1_L001.fastq.gz`` and ``S10_001_R2_L001.fastq.gz`` in the folder defined by ``link_directory`` folder will be matched to the sample name ``S10``.

If needed, an `OldSampleName` column can be added to the file, when the read file names and the desired new sample names can not be match simply by testing the identity at the start of both strings. 

.. csv-table:: Local data example with old sample names
   :header: "SampleName", "ScientificName", "OldSampleName"
   
   "S10","Staphylococcus aureus","Staaur-10"
   "S1","Staphylococcus aureus","Staaur-1"	
   "S2","Staphylococcus aureus","Staaur-2"	
   "S3","Staphylococcus aureus","Staaur-3"	
   "S4","Staphylococcus aureus","Staaur-4"	
   "S5","Staphylococcus aureus","Staaur-5"	
   "S6","Staphylococcus aureus","Staaur-6"	
   "S7","Staphylococcus aureus","Staaur-7"	
   "S8","Staphylococcus aureus","Staaur-8"	
   "S9","Staphylococcus aureus","Staaur-9"
      
In this case, the files ``Staaur-10_S10_L001_R1_001.fastq.gz`` and ``Staaur-10_S10_L001_R2_001.fastq.gz`` in the folder defined in ``link_directory`` will be matched to the sample name ``S10``.

   
SRA samples
-----------

SRA samples will be determined based on the tabulated file whose full path must be passed to the variable ``sra_samples``. The RunInfo files that can be downloaded through the NCBI can be directly passed without any modification. Otherwise, four columns must be defined.

.. csv-table:: SRA data example
   :header: "Run","SampleName", "LibraryLayout", "ScientificName"
	 
   "ERR1140788","Mycobacterium_tuberculosis_N0145-Lineage_2","paired","Mycobacterium tuberculosis"
   "SRR006916","Mycobacterium_tuberculosis_K21-Lineage_1","single","Mycobacterium tuberculosis"
   "SRR022880","Mycobacterium_tuberculosis_africanum_4141_04-Lineage_6","single","Mycobacterium tuberculosis"
   "ERR551336","Mycobacterium_tuberculosis_11821-03-Lineage_5","paired","Mycobacterium tuberculosis"
   "ERR756345","Mycobacterium_tuberculosis_aethiop_vetus-Lineage_7","paired","Mycobacterium tuberculosis"
   "ERR238745","Mycobacterium_tuberculosis_SG1-Lineage_3","single","Mycobacterium tuberculosis"
