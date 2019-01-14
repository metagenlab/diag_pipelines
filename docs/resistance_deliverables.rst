
------------------
rgi and mykrobe analysis
------------------

* ``report/resistance/rgi_summary.xlsx``: summary of RGI results for every samples (one per sheet)
* ``report/resistance/mykrobe_summary.xlsx``: summary of Mykrobe results for every samples (one per sheet, only for `Mycobacterium` genus or `Staphylococcus aureus`)
* ``report/resistance/{sample}_rgi_report.html``: summary of RGI results as html document with cross references to the CARD database

------------------
*Mycobacterium tuberculosis* SNPs variant associated to resistance
------------------

* ``sample/{sample_name}/resistance/bwa/miotto_high_moderate_minimum_confidence_annotated/mutations.vcf``: VCF files of all markers in `Miotto et al. 2017 European Respiratory Journal <http://erj.ersjournals.com/content/50/6/1701354>`_
* ``sample/{sample_name}/resistance/bwa/mykrobe_annotated/mutations.vcf``: VCF files of all markers in `Bradley et al. 2015 Nature Communications <https://www.nature.com/articles/ncomms10063>`_
* ``sample/{sample_name}/resistance/bwa/walker_resistant_annotated/mutations.vcf``: VCF files of all markers in `Walker et al. 2015 Lancet Infectious Diseases <https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(15)00062-6/abstract>`_

* ``sample/{sample_name}/resistance/bwa/rgi_annotated_full_2_0_0/mutations.vcf``: VCF files of all markers in version 2.0.0 of `the CARD database <https://card.mcmaster.ca/>`_
