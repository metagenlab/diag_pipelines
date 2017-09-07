# Requirements
Miniconda  
Snakemake > 4.0.0  
Setting kraken database  
```
        wget -O - https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz > krakendb.tgz
        tar xzf krakendb.tgz
        mv minikraken_* krakendb/
        rm -rf krakendb.tgz
```
Illumina adapter file for trimming  
Modify config file to your needs
The paired reads files have to be uncompressed, stored in the folder reads/raw/ from where the snakemake command is called and be their names have be ${ISOLATE_NAME}_R1_001.fastq and ${ISOLATE_NAME}_R2_001.fastq.

Call Snakemake like this:
```
snakemake -s ${PATH_TO_ASSEMBLY_RULES} --use-conda --conda-prefix ${PATH_WHERE_CONDA_ENVS_WILL_BE_STORED} ${ISOLATE_NAME}_report.html
```