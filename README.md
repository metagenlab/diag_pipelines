# Pipeline for assembling a set of pair read sequecing files and performing comparative genomics analysis

This pipeline is composed of two modules. The first one assembles paired reads, checks the identity of the species, annotates the genome, controls for contamination and generates an EMBL file ready to be submitted to the European Nucleotide Archive.

The second module uses the assembled genomes, and downloads from requested species or genera already sequenced genomes to perform an orthologous gene search, an ANI calculation, and a phylogeny based on one to one orthologous protein coding genes. It accomodates any number of sequenced genome. 

The pipeline uses the conda package manager, and will run provided snakemake, Miniconda3 are available on the osX or Linux computer.
Setting up a minikraken database is also needed : 
```
wget -O - https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz > krakendb.tgz
tar xzf krakendb.tgz
mv minikraken_* krakendb/
rm -rf krakendb.tgz
```

The paired reads files have to be uncompressed, stored in the folder reads/raw/ relative from where the snakemake command is called and their names be `${ISOLATE_NAME}_R1_001.fastq` and `${ISOLATE_NAME}_R2_001.fastq`. For each isolate, a file called ${ISOLATE_NAME}.yaml must be present in the main folder (see example in `strain_name.yaml`)

Modify config.yaml and copy it in the folder from which the analysis will be run. For the RefSeq genome search, do not include large genera or too many genomes will be added.   

Check every thing is order by performing first a dry run :

```
snakemake --snakefile ${PATH_TO_COMPARATIVE_RULES_FILE} --dryrun --use-conda --conda-prefix=${PATH_WHERE_CONDA_ENV_WILL_BE_STORED} report.html
```

A directed acyclic graph (dag) of the run can be generated to visualize all jobs that will be performed :

```
snakemake --snakefile ${PATH_TO_COMPARATIVE_RULES_FILE} --use-conda --conda-prefix=${PATH_WHERE_CONDA_ENV_WILL_BE_STORED} --dag report.html | dot -Tsvg > dag.svg
```

Finally run the full pipeline : 

```
snakemake --snakefile ${PATH_TO_COMPARATIVE_RULES_FILE} --use-conda --conda-prefix=${PATH_WHERE_CONDA_ENV_WILL_BE_STORED} report.html 
```

You can find in `dag.svg` an example of the pipeline with 4 strains that are assembled and then compared to RefSeq genomes.
