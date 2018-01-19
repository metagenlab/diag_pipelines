# DEPENDENCIES
  Snakemake (http://snakemake.readthedocs.io/en/latest/) 4.0.0 and up
  
  Miniconda (https://conda.io/miniconda.html) preferably version with python 3.6



# GENERAL USE

```
snakemake --snakefile ${PATH_TO_GIT_FOLDER}/workflows/typing/general_workflow.rules --use-conda --conda-prefix ${PATH_TO_CONDA_INSTALLATION} --configfile config.yaml
```


Update the config file for your needs.

# GENERATING FILES OF INTEREST

The pipeline works by asking the generation of the files of interest for a particular analysis. This assumes knowing what their locations is going to be.

Some examples:

```
snakemake --snakefile ${PATH_TO_GIT_FOLDER}/workflows/typing/general_workflow.rules --use-conda --conda-prefix ${PATH_TO_CONDA_INSTALLATION} --configfile config.yaml quality/multiqc_report.html
```

This will assemble and annotate every samples present in the `links` folder, and generate a multiqc report.


```
snakemake --snakefile ${PATH_TO_GIT_FOLDER}/workflows/typing/general_workflow.rules --use-conda --conda-prefix ${PATH_TO_CONDA_INSTALLATION} --configfile config.yaml virulence_summary.xlsx
```

This will generate a summary excel file for the virulence factors of the strains in the `links` folder, based on the virulence factors annotated in the file defined on the config file.



```
snakemake --snakefile ${PATH_TO_GIT_FOLDER}/workflows/typing/general_workflow.rules --use-conda --conda-prefix ${PATH_TO_CONDA_INSTALLATION} --configfile config.yaml resistance_summary.xlsx
```

This will generate a summary excel file for the resistance factors of the strains in the `links` folder, using the softwares defined in the config file.
