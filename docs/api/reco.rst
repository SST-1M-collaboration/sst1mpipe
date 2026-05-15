Reconstruction
==============

Overview
--------

The :mod:`sst1mpipe.reco` namespace contains the machine-learning and geometric
reconstruction logic used to transform calibrated event parameters into DL2
quantities. It is responsible for training the random-forest models used by the
pipeline and for combining mono or stereo telescope information into final
event-level estimates.

Key responsibilities
--------------------

* train DISP, energy, particle-classification, and misdirection models;
* apply trained models to DL1 tables and derive reconstructed event features;
* combine telescope-level predictions into stereo event quantities.

Usage example
-------------

Training is typically done on prepared gamma and proton MC samples:

.. code-block:: python

    from sst1mpipe.reco import train_models

    train_models(
        params_gamma=gamma_training_sample,
        params_protons=proton_training_sample,
        config=config,
        outdir="models",
        telescope="tel_001",
        plot=True,
    )

For stereo reconstruction, the reconstructed telescope rows are merged back into
an event-level table:

.. code-block:: python

    from sst1mpipe.reco import stereo_reconstruction

    dl2_stereo = stereo_reconstruction(
        params=stereo_params,
        config=config,
        ismc=False,
        telescopes=["tel_021", "tel_022"],
    )

Architecture notes
------------------

This module sits between feature engineering and scientific validation. It
expects feature-rich DL1 tables, applies scikit-learn models, and returns
tables that are later consumed by :mod:`sst1mpipe.performance` and
:mod:`sst1mpipe.analysis`. Training can be both CPU- and memory-intensive,
especially when temporary models are created to enrich classifier features.

API reference
-------------

.. automodule:: sst1mpipe.reco
   :members:
   :undoc-members:
   :show-inheritance:
