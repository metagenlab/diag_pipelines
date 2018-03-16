.. _resistance:

==========
Resistance
==========

Depends on the :ref:`assembly_quality` workflow.

----------
Parameters
----------

* ``resistance_prediction_softwares``: list of software for genetic resistance assessment. Possible values: ``mykrobe`` and ``rgi``.
  
.. * ``currated_resistance_genes``: file of trusted genes involved in resistance. An example is available in the folder ``data/mycobacterium/db/``

------------
Deliverables
------------

* ``samples/{sample_name}/annotation/resistance/rgi.tsv``: results files for RGI 
* ``samples/{sample_name}/annotation/resistance/mykrobe.tsv``: results file for mykrobe
