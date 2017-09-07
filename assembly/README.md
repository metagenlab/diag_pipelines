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
Modify config file  