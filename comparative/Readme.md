# Comparative analysis with genomes from RefSeq
Modify comparative_search.yaml to search RefSeq for species/genus. Do not include large genera or too many genomes will be added.
Include genomes and proteomes to be added in local_data/genomes and local_data/proteomes/to_be_added/. The extensions of the files have to be .fna and .faa.
Include a list of the ids of the isolates in local_data/ids.txt.   
Call snakemake:   
```
snakemake -s ${PATH_TO_COMPARATIVE_RULES} --use-conda --conda-prefix=${PATH_WHERE_CONDA_ENV_WILL_BE_STORED} report.html
```