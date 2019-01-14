

~~~~~~~~~~~~~~~~~~~~~~
Assembly quality
~~~~~~~~~~~~~~~~~~~~~~

* ``report/multiqc_assembly/multiqc_report.html``: quality control report based on the results of **fastqc**, **trimmomatic**, **qualimap**, **quast** and **prokka** for every sample
* ``report/contamination/low_coverage_contigs/{sample}.html``: quality control report based on the results of **fastqc**, **trimmomatic**, **qualimap**, **quast** and **prokka** for every sample
* ``report/multiqc_assembly/multiqc_report.html``: quality control report based on the results of **fastqc**, **trimmomatic**, **qualimap**, **quast** and **prokka** for every sample

~~~~~~~~~~~~~~~~~~~~~~
Assembly and annotation
~~~~~~~~~~~~~~~~~~~~~~

* ``samples/{sample}/assembly/spades/``: folder containing all files from **SPADES** assembly (contigs.fasta, assembly_graph.fastg, contigs.paths,...)
* ``samples/{sample_name}/annotation/``: folder containing all annotation files from **prokka**
