# DEPENDENCIES
  Docker  

# GENERAL USE
Once you have pulled the docker image on your computer: 
```
docker run --rm --mount source="$(pwd)",target=/home/pipeline_user,type=bind diag_pipeline sh -c 'snakemake --snakefile $pipeline_folder/workflows/general_workflow.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml'
```

Update the config file for your needs.

# GENERATING FILES OF INTEREST

The pipeline works by asking the generation of the files of interest for a particular analysis. This assumes knowing what their locations is going to be.


Some examples:

```
docker run --rm --mount source="$(pwd)",target=/home/pipeline_user,type=bind diag_pipeline sh -c 'snakemake --snakefile $pipeline_folder/workflows/general_workflow.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml quality/multiqc/self_genome/multiqc_report.html'
```

This will assemble and annotate every samples present in the `links` folder, and generate a multiqc report.


```
docker run --rm --mount source="$(pwd)",target=/home/pipeline_user,type=bind diag_pipeline sh -c 'snakemake --snakefile $pipeline_folder/workflows/general_workflow.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml virulence_summary.xlsx'
```

This will generate a summary excel file for the virulence factors of the strains in the `links` folder, based on the virulence factors annotated in the file defined on the config file.



```
docker run --rm --mount source="$(pwd)",target=/home/pipeline_user,type=bind diag_pipeline sh -c 'snakemake --snakefile $pipeline_folder/workflows/general_workflow.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml resistance_summary.xlsx'
```

This will generate a summary excel file for the resistance factors of the strains in the `links` folder, using the softwares defined in the config file.


```
docker run --rm --mount source="$(pwd)",target=/home/pipeline_user,type=bind diag_pipeline sh -c 'snakemake --snakefile $pipeline_folder/workflows/general_workflow.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml typing/freebayes/core_parsnp/34528/bwa/distances_from_merged_pairs_of_vcf.xlsx'
```

This will generate a snp-distance matrix of all samples present in `links`, only on the core genome calculated with parsnp and with all complete genomes of the species defined in the `taxid` variable of the config file, mapped on the assembly (from https:/www.ncbi.nlm.nih.gov/assembly/) whose `id` is 34528 (nctc 8325, *Staphylococcus aureus* reference genome) with bwa


```
docker run --rm --mount source="$(pwd)",target=/home/pipeline_user,type=bind diag_pipeline sh -c 'snakemake --snakefile $pipeline_folder/workflows/general_workflow.rules --use-conda --conda-prefix $conda_folder --configfile config.yaml typing/mlst/summary.xlsx'
```

This will generate an Excel summary file of the MLST of all samples present in `links`, based on the software mlst (https:/github.com/tseemann/mlst)

