Calibration
===========

Overview
--------

The :mod:`sst1mpipe.calib` namespace implements calibration steps that adjust
raw SST-1M camera data before higher-level reconstruction. These routines are
used close to the R0-to-R1 and early DL1 stages, where waveform corrections and
instrument response assumptions have the largest impact on downstream quality.

Key responsibilities
--------------------

* load or derive window-transmittance corrections;
* apply PDE- and saturation-related calibration factors;
* provide the calibrator class used by the lower-level processing scripts.

Usage example
-------------

.. code-block:: python

    from sst1mpipe.calib import get_window_corr_factors

    window_corr, window_file = get_window_corr_factors(
        telescope=21,
        config=config,
    )

Architecture notes
------------------

Calibration is deliberately isolated from the later machine-learning and
science-analysis layers. The functions in this module generally consume config
objects and array data, then return corrected values without orchestrating the
rest of the pipeline. This keeps low-level telescope-specific adjustments
separate from the higher-level event classification logic.

API reference
-------------

.. automodule:: sst1mpipe.calib
   :members:
   :undoc-members:
   :show-inheritance:
