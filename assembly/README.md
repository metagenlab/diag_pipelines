# Requirements
Modify strain config file (.yaml) to your needs and place it in the folder where snakemake is going to be called.

Call Snakemake like this:
```
snakemake -s ${PATH_TO_ASSEMBLY_RULES_FILE} --use-conda --conda-prefix ${PATH_WHERE_CONDA_ENVS_WILL_BE_STORED} ${ISOLATE_NAME}/report.html
```
