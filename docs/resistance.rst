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


--------------------------------------------
Mycobacterium tuberculosis specific analyses
--------------------------------------------

If ``species`` is defined to Mycobacterium_tuberculosis, specific variant markers of resistance can be searched by mapping to the genome of H37Rv, genotyping with GATK and annotation of the resulting VCF. Different databases of markers can be searched:
* Markers from `Walker et al. 2015 Lancet Infectious Diseases <https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(15)00062-6/abstract>`_
* Markers from `the CARD database <https://card.mcmaster.ca/>`_
* Markers from `Miotto et al. 2017 European Respiratory Journal <http://erj.ersjournals.com/content/50/6/1701354>`_
* Markers derived `Bradley et al. 2015 Nature Communications <https://www.nature.com/articles/ncomms10063>`_
 

------------
Deliverables
------------

.. include:: resistance_deliverables.rst
