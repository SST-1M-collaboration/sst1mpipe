Analysis
========

Overview
--------

The :mod:`sst1mpipe.analysis` namespace contains the post-reconstruction tools
that help transform DL2 and DL3 event tables into analysis-ready coordinate
representations, source-centric features, and publication-facing plots. This is
the layer that connects raw reconstructed quantities to higher-level scientific
interpretation.

Key responsibilities
--------------------

* convert reconstructed alt/az coordinates into sky coordinates;
* derive source-relative quantities such as wobble labels and camera-frame
  positions;
* provide plotting helpers for count maps, angular resolution, sensitivity, and
  other analysis products.

Usage example
-------------

The functions below are typically used after DL2 reconstruction, once the event
table already contains reconstructed directions and telescope pointing
information:

.. code-block:: python

    from sst1mpipe.analysis import add_reco_ra_dec
    from sst1mpipe.utils import get_horizon_frame

    horizon_frame = get_horizon_frame(
        config=config,
        telescope="tel_021",
        times=dl2_data["local_time"],
    )
    dl2_with_radec = add_reco_ra_dec(dl2_data, horizon_frame=horizon_frame)

Architecture notes
------------------

These helpers sit close to the scientific analysis boundary. They do not train
models or read telescope files directly; instead, they operate on already
prepared event tables and coordinate frames. Most functions copy or enrich the
input tables, so memory usage can become noticeable when large DL2 or DL3
samples are processed interactively in notebooks.

API reference
-------------

.. automodule:: sst1mpipe.analysis
   :members:
   :undoc-members:
   :show-inheritance:
