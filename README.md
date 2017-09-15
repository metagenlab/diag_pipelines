# Pipeline for assembling a set of pair read sequecing files and performing comparative genomics analysis.

This pipeline is composed of two modules. The first one assembles paired reads, checks the identity of the species, annotates the genome, controls for contamination and generates an EMBL file ready to be submitted to the European Nucleotide Archive.

The second pipeline uses the assembled genomes, and downloads from requested species or genera already sequenced genomes to perform an orthologous gene search, an ANI calculation, and a phylogeny based on 1 to 1 orthologous protein coding genes. It accomodates any number of sequenced genome, modification of comparative/config.yaml and config file for each isolate only is needed.

It uses the conda package manager, and will run provided Snakemake, Miniconda3 are available on the osX or Linux computer.
Setting up a minikraken database is also needed : 
```
        wget -O - https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz > krakendb.tgz
        tar xzf krakendb.tgz
        mv minikraken_* krakendb/
        rm -rf krakendb.tgz
```

The paired reads files have to be uncompressed, stored in the folder reads/raw/ relative from where the snakemake command is called and their names be ${ISOLATE_NAME}_R1_001.fastq and ${ISOLATE_NAME}_R2_001.fastq.  

Modify config.yaml to search RefSeq for species/genus. Do not include large genera or too many genomes will be added. 
