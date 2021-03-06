.. _mykrobe:

=======
Mykrobe
=======

`Mykrobe <http://www.mykrobe.com/products/predictor/>`_ can be used on *Staphylococcus aureus* and *Mycobacterium tuberculosis* samples only. Predictions are based directly on the fastq reads. Version ``0.5.6`` is used.

----------
Parameters
----------
For *M. tuberculosis*, two different panels of mutations can be analysed, by defining the variable ``mykrobe_panel``. Two different values are possible:

- ``bradley-2015``, from  `Bradley et al. 2015, Nature Communications <http://www.mykrobe.com/wp-content/uploads/2014/04/ncomms10063.pdf>`_ (default value)
- ``walker-2015``, from `Walker et al. 2015, Lancet Infectious Diseases <https://www.ncbi.nlm.nih.gov/pubmed/26116186>`_

The threshold of confidence for Mykrobe can also be defined with ``mykrobe_confidence``. The default value is ``10``.
