Utilities
=========

Overview
--------

The :mod:`sst1mpipe.utils` namespace groups reusable helper functions that are
shared across the pipeline. These functions handle coordinate conversions,
feature engineering, telescope metadata access, image cleaning setup, and
various domain-specific convenience operations used by scripts and notebooks.

Key responsibilities
--------------------

* provide shared event- and table-level helpers used by multiple modules;
* centralize common telescope and coordinate utilities;
* expose image-cleaning and NSB-related support functions.

Usage example
-------------

.. code-block:: python

    from sst1mpipe.utils import get_telescopes, get_wr_timestamp

    telescopes = get_telescopes(config=config)
    timestamp = get_wr_timestamp(event)

Architecture notes
------------------

This namespace intentionally collects code that does not belong to a single
pipeline stage but would otherwise be duplicated. Because these helpers are used
widely, changes here can have broad effects across calibration, reconstruction,
and analysis workflows. Prefer keeping interfaces stable and documenting
cross-cutting assumptions clearly.

API reference
-------------

.. automodule:: sst1mpipe.utils
   :members:
   :undoc-members:
   :show-inheritance:
