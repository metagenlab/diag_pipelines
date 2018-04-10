.. _rgi:


===================================================
Comprehensive Antibiotic Resistance Database (CARD)
===================================================


The Resistance Gene Identifier (RGI) from the `Comprehensive Antibiotic Resistance Database <https://card.mcmaster.ca/home>`_, version ``3.2.1`` is used, with data version ``1.1.9``. Predictions are based on the assembled contigs.


Ontology
========

The CARD ontology is parsed after the analysis are run to summarize the antibiotics for which a gene confering resistance has been identified (testing for the relationship equals to ``confers_resistance_to_drug`` or ``confers_resistance_to``). Find the result of the parsing for each sample in:

* ``samples/{sample}/resistance/rgi_ontology.xlsx``
