Instrument
==========

Overview
--------

The instrument modules describe the SST-1M camera layout and the geometry
objects needed by calibration and visualization code. They bridge external
camera descriptions from the CTS tooling to the geometry abstractions used by
``ctapipe`` and the rest of this repository.

Key responsibilities
--------------------

* wrap the CTS camera description with SST-1M-specific defaults;
* generate pixel geometry and adjacency structures;
* expose reusable camera metadata to higher-level modules.

Usage example
-------------

.. code-block:: python

    from sst1mpipe.instrument.camera import Camera

    camera = Camera()
    geometry = camera.geometry

Architecture notes
------------------

These modules sit on an integration boundary between SST-1M-specific resources
and external packages such as ``cts_core`` and ``ctapipe``. The geometry
helpers are read-mostly and are typically instantiated once per workflow, so
clarity and correctness matter more than micro-optimizing their runtime.

API reference
-------------

.. automodule:: sst1mpipe.instrument

Camera helpers
--------------

.. automodule:: sst1mpipe.instrument.camera
   :members:
   :undoc-members:
   :show-inheritance:

Geometry helpers
----------------

.. automodule:: sst1mpipe.instrument.geometry
   :members:
   :undoc-members:
   :show-inheritance:
