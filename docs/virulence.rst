.. _virulence:
  
Virulence
=========

Depends on the :ref:`assembly_quality` workflow.

----------
Parameters
----------

* ``virulence_factors``: file with list of uniprot accession of virulence factors. An example is available in the folder ``data/staph/db/``
* ``virulence_percentage_identity_cutoff``: amino acid identity cut off for considering a match
* ``virulence_coverage_cutoff``: coverage cut off for considering a match

------------
Deliverables
------------

*  ``virulence_summary.xlsx``: summary of virulence proteins found in every samples.
