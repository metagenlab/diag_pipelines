## Pipeline benchmark tools


# General Information

The python script present in this folder is meant to validate pipelines with a set of reference data in addition to comparing the output of different pipelines between them.
The goal is to ensure that new editions of the diag_pipelines tool used for routine sequencing of clinically relevant strains are benchmarked reproducibly on the same reference dataset.

This folder contains the following files:

The two tables containing the raw data for the reference data for both gram negative (reference_table_gramneg.tsv) and gram positive samples (reference_table_grampos.tsv). Gram negative samples contain 
MIC values in addition to microarray data and were taken from the following study: https://pubmed.ncbi.nlm.nih.gov/24930405/ (the Vogne_raw_data.xlsx file contains the raw data as provided by the authors).
Gram positive samples contain PCR results for 2 targeted resistance genes.

The samples folder contains the config file and tsv file with the name and location of the samples. Those need to be used to run the pipeline and obtain the tsv files necessary for the analysis.

The examples folder contains the files outputed by running the script using the combined_detail_NR.tsv file

The diag_pipelines_2.6 contains the previous iteration of the pipeline as is used by default if no other tsv was provided.
 
The code outputs two files: 

A pdf file with the specificity and sensitivity of two previous pipeline iterations in addition to the latest pipeline compared to reference MIC phenotype, PCR and Microarray data.

An html table that shows the differences in prediction between the previous iteration of the pipeline and the new one. The first column corresponds to the genes predicted
in the older iteration and not in the newer one and the second column shows the genes predicted in the new iteration and not the older one.

Dataset used:

A set of gram negative samples from the following study (https://pubmed.ncbi.nlm.nih.gov/24930405/)  were used used for which there was both MIC Phenotype, Microarray and Whole Genome Sequencing data.
For gram positive samples, a set of internally produced Whole Genome Sequencing data with PCR confirmations of resistance genes was used made up 20 samples that were sequenced twice.

Phenotype Sensitivity and Specificity score:

True Positive: If a sample had a resistance MIC for at least one carbapenem as per eucast Breakpoint table and at least one carbapenemase gene was detected
False Positive: If a sample had no MIC associated with carbapenem resistant and at least one carbapenemase gene was detected
True Negative: If a sample had no MIC associated with carbapenem resistant and no carbapenemase gene was detected
False Negative: If a sample had a resistance MIC for at least one carbapenem as per eucast Breakpoint table and no carbapenemase gene was detected

Microarray Sensitivity:

Considering that Whole Genome Sequencing has greater discriminatory power than the targeting of individual genes by Microarray, the goal was to verify that the gene families detected
by microarray are also preticted by the pipeline.

True Positive: If a gene is detected by Microarray and also predicted in the pipeline output
False Negative: If a gene is detected by Microarray and not predicted in the pipeline output

PCR Sensitivity and Specificity score:

True Positive: If a sample had a confirmed PCR assay for one of the two targeted genes and that same gene was predicted in the pipeline output.
False Positive: If a sample had no confirmed PCR assay for a targeted gene but a gene was predicted in the pipeline output.
True Negative: If a sample had no confirmed PCR assay for one of the two targeted genes and no gene was predicted in the pipeline output.
False Negative: If a sample had a confirmed PCR assay for one of the two targeted genes and no gene was predicted in the pipeline output.

# Installation:

Creating conda env:
```
conda create --name pipeline_benchmark python=3.9
conda activate pipeline_benchmark
pip install -r requirements.txt
```

Install wkhtmltopdf:

```
apt-get install wkhtmltopdf
```

Usage:

To run this script, you need a tsv file with the gene predictions as outputed by the diag_pipelines.

```
python benchmark_calculations.py -h
```
Example of commands

```
python benchmark_calculations.py combined_detail_NR.tsv -p diag_pipelines_2.6.0/combined_detail_NR_2.6.tsv -o results
```

# Output

pipeline_benchmark_summary.pdf:

a PDF file with tables containing different statistics of the previous and current pipeline iteration

pipeline_comparison.html:

an HTML file with a prediction comparison of the two pipelines





