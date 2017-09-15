# Comparative analysis with genomes from RefSeq, assembles the local genomes 

Modify config.yaml to search RefSeq for species/genus. Do not include large genera or too many genomes will be added. 
Include all the name of the strains for which you have sequencing read files. 

Call snakemake:

```
snakemake -s ${PATH_TO_COMPARATIVE_RULES} --use-conda --conda-prefix=${PATH_WHERE_CONDA_ENV_WILL_BE_STORED} report.html
```
