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

* ``samples/{sample_name}/annotation/resistance/rgi.xlsx``: results files for RGI, for each sample
* ``samples/{sample_name}/annotation/resistance/mykrobe.xlsx``: results file for mykrobe, for each sample
* ``resistance/{mykrobe_or_rgi}_summary.xlsx``: summary for every samples (one per sheet)  
