sst1mpipe package
=================

Overview
--------

``sst1mpipe`` is the Python implementation of the SST-1M prototype processing
pipeline. It spans low-level calibration, reconstruction, scientific
performance studies, and production utilities that convert raw or simulated
events into higher-level data products suitable for downstream analysis.

Key responsibilities
--------------------

* define the public package layout used by scripts and notebooks;
* group the pipeline into calibration, I/O, reconstruction, and performance
  stages;
* centralize telescope-specific helpers that are reused across production and
  exploratory workflows.

Typical workflow
----------------

A common developer workflow starts by loading a JSON configuration and DL1 or
DL2 tables from :mod:`sst1mpipe.io`, performing reconstruction or feature
selection through :mod:`sst1mpipe.reco` and :mod:`sst1mpipe.analysis`, and then
computing validation products with :mod:`sst1mpipe.performance`.

Subpackages
-----------

The package is organized around the main pipeline stages:

* :mod:`sst1mpipe.calib` for calibration and waveform corrections;
* :mod:`sst1mpipe.io` for reading and writing pipeline products;
* :mod:`sst1mpipe.reco` for event reconstruction and model training;
* :mod:`sst1mpipe.performance` for performance studies and sensitivity;
* :mod:`sst1mpipe.analysis` for higher-level analysis helpers and plots;
* :mod:`sst1mpipe.utils` for shared utilities used across the codebase;
* :mod:`sst1mpipe.instrument` for camera and geometry descriptions.

API reference
-------------

.. automodule:: sst1mpipe
