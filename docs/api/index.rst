API reference
=============

Overview
--------

This section gathers the public Python interfaces that make up the
``sst1mpipe`` analysis chain. The reference is intended to complement the
workflow-oriented narrative pages by explaining where each module fits in the
pipeline and by pointing to the most relevant entry points for developers who
need to extend, debug, or automate the project.

How to use this reference
-------------------------

Start with :doc:`package` for a high-level map of the package, then continue
with the module pages that match the pipeline stage you are working on:

* :doc:`calib` for raw waveform and calibration corrections.
* :doc:`io` for configuration handling and DL1/DL2/DL3 table access.
* :doc:`reco` for machine-learning training and mono/stereo reconstruction.
* :doc:`performance` for sensitivity curves and validation products.
* :doc:`analysis` for higher-level event selection and plotting helpers.
* :doc:`instrument` and :doc:`utils` for shared support code.

Architecture notes
------------------

The package is organized around data-level transitions. Most user-facing
scripts combine these namespaces rather than reimplementing the pipeline logic:
configuration is loaded through :mod:`sst1mpipe.io`, event features are
produced and transformed in :mod:`sst1mpipe.reco` and :mod:`sst1mpipe.analysis`,
and validation products are computed in :mod:`sst1mpipe.performance`.

.. toctree::
   :maxdepth: 1

   package
   analysis
   calib
   instrument
   io
   performance
   reco
   utils
