.. _resistance:

==========
Resistance
==========

Depends on the :ref:`epidemiology` workflow (for genotyping samples)

----------
Parameters
----------

* ``resistance_prediction_softwares``: list of software for genetic resistance assessment. Possible values: ``mykrobe`` and ``rgi``.
  
-------------------
Available softwares
-------------------

.. toctree::
   resistance_softwares/rgi.rst
   resistance_softwares/mykrobe.rst



------------
Deliverables
------------

* ``samples/{sample_name}/annotation/resistance/rgi.tsv``: results files for RGI 
* ``samples/{sample_name}/annotation/resistance/mykrobe.tsv``: results file for mykrobe
