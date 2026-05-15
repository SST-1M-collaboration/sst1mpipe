Performance
===========

Overview
--------

The :mod:`sst1mpipe.performance` namespace evaluates the quality of the
reconstruction chain. It computes quantities such as energy resolution, angular
resolution, ROC curves, and flux sensitivity, and it stores the intermediate
cuts needed by later IRF and DL3 production steps.

Key responsibilities
--------------------

* validate trained models on gamma and proton test samples;
* derive sensitivity curves and energy-dependent event-selection cuts;
* generate performance plots and HDF tables for later pipeline stages.

Usage example
-------------

The module is usually driven by testing DL2 samples for simulated gammas and
protons:

.. code-block:: python

    from sst1mpipe.performance import evaluate_performance, sensitivity

    evaluate_performance(
        gamma_file="gamma_test_dl2.h5",
        proton_file="proton_test_dl2.h5",
        outdir="performance",
        config=config,
        telescope="tel_001",
        save_fig=True,
    )

    sensitivity(
        "gamma_test_dl2.h5",
        "proton_test_dl2.h5",
        outdir="performance",
        config=config,
        telescope="tel_001",
        gammaness_cuts="significance",
        theta2_cuts="efficiency",
    )

Architecture notes
------------------

These routines sit at the scientific validation boundary of the pipeline. They
read reconstructed event tables, apply weighting assumptions from the configured
spectra, and generate summary products that are reused in offline analyses.
Most computations scale with the size of the DL2 samples and can be expensive
when repeated interactively, especially when energy-dependent cuts and multiple
plots are requested.

API reference
-------------

.. automodule:: sst1mpipe.performance
   :members:
   :undoc-members:
   :show-inheritance:
