Input and output
================

Overview
--------

The :mod:`sst1mpipe.io` namespace defines the repository's data access layer.
It loads configuration files, reads SST-1M DL1 and DL2 tables, merges auxiliary
tables, and writes derived products back to HDF5. Most command-line scripts and
notebooks enter the pipeline through this module.

Key responsibilities
--------------------

* resolve default data and Monte Carlo configuration files;
* read DL1, DL2, monitoring, and photon-list products into Astropy or pandas
  tables;
* write derived products such as DL2 event tables and supporting metadata;
* keep table layout and path conventions consistent across the pipeline.

Usage example
-------------

The snippet below shows a typical mono-analysis setup that loads the default
configuration, reads a DL1 table, and writes a derived DL2 product:

.. code-block:: python

    from sst1mpipe.io import load_config, load_dl1_sst1m, write_dl2

    config = load_config(None, ismc=True)
    dl1 = load_dl1_sst1m(
        "gamma_test_dl1.h5",
        tel="tel_001",
        config=config,
        table="pandas",
        quality_cuts=True,
    )
    write_dl2(dl1, output_file="gamma_test_dl2.h5", telescope="tel_001", config=config)

Architecture notes
------------------

This namespace forms the boundary between on-disk SST-1M products and the
in-memory tables used by reconstruction and analysis modules. Many functions
return either Astropy or pandas objects depending on the workflow, so callers
should choose the representation that matches the next processing step.
Converting large tables to pandas improves interoperability with scikit-learn
but increases memory pressure and removes some multi-dimensional columns.

API reference
-------------

.. automodule:: sst1mpipe.io
   :members:
   :undoc-members:
   :show-inheritance:
