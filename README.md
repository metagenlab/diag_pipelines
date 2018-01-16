# GENERAL USE

```
snakemake --snakefile ${PATH_TO_GIT_FOLDER}/workflows/typing/general_workflow.rules --use-conda --conda-prefix ${PATH_TO_CONDA_INSTALLATION} --configfile config.yaml
```


Update the config file for your needs.

In the `general_workflow.rules` file, change the value of the snakemake_path variable to the full path of the git folder.
