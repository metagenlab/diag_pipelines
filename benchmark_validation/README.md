## Pipeline benchmark tools

The python present in this folder is meant to validate pipelines with a set of reference data in addition to comparing the output of different pipelines between them.
The goal is to ensure that new editions of the diag_pipelines tool used for routine sequencing of clinically relevant strains are benchmarked reproducibly on the same reference dataset.

The code outputs two files: 

A pdf file with the specificity and sensitivity of two previous pipeline iterations in addition to the latest pipeline compared to reference MIC phenotype data.

An html table that shows the differences in prediction between the previous iteration of the pipeline and the new one. The first column corresponds to the genes predicted
in the older iteration and not in the newer one and the second column shows the genes predicted in the new iteration and not the older one.

Phenotype Sensitivity and Specificity score:

True Positive: If a sample had a resistance MIC for at least one carbapenem as per eucast Breakpoint table and at least one carbapenemase gene was detected
False Positive: If a sample had no MIC associated with carbapenem resistant and at least one carbapenemase gene was detected
True Negative: If a sample had no MIC associated with carbapenem resistant and no carbapenemase gene was detected
False Negative: If a sample had a resistance MIC for at least one carbapenem as per eucast Breakpoint table and no carbapenemase gene was detected

Microarray Sensitivity:


