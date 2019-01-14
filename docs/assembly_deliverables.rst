* ``report/multiqc_assembly/multiqc_report.html``: quality control report based on the results of **fastqc**, **trimmomatic**, **qualimap**, **quast** and **prokka** for every sample
* ``samples/{sample_name}/annotation/``: folder containing all annotation files from **prokka**
* ``report/contamination/mash/assembly/distances_formated.xlsx``: mash distances to all bacterial, viral and human RefSeq genomes, ordered by significance, for each sample (mash of the assembly)
* ``report/contamination/mash/reads/distances_formated.xlsx``: mash distances to all bacterial, viral and human RefSeq genomes, ordered by significance, for each sample (mash of the reads)
* ``report/contamination/centrifuge/{sample_name}/formatted.tsv``: summary of the taxonomic classification of trimmed reads with **centrifuge**
