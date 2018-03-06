.. diag_pipelines documentation master file, created by
   sphinx-quickstart on Wed Feb 28 12:10:42 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========================================================
Documentation for the genomics and metagenomics workflows
=========================================================

Routine procedures for diagnostic purposes using microbial genomics and metagenomics.

Workflows for epidemiology, anti-microbial resistance genotyping and virulence factors identification have been implemented using the `Snakemake <http://snakemake.readthedocs.io/en/stable/>`_ workflow management system with `bioconda <https://bioconda.github.io/>`_ integration for software dependency. `Docker <https://hub.docker.com/r/metagenlab/diag_pipelines/>`_ images of main releases are available.

As a general rules, any ``variable`` referenced in this documentation must be either:

* Defined in the yaml config file that is passed to snakemake by ``--configfile``
* Defined directly in the snakemake command by ``--config variable=$value``

.. toctree::
   general_processes.rst
   workflows.rst
   :maxdepth: 2
   :caption: Contents:
   



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
