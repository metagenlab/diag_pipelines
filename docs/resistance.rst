.. _resistance:

==========
Resistance
==========

Depends on the :ref:`epidemiology` workflow (for genotyping samples)

----------
Parameters
----------

* ``species``: Species name to determine which software are available to run for your sample


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
