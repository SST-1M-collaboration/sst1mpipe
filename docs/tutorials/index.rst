Tutorials
=========

Examples of how to use sst1mpipe, as Jupyter notebooks. Each notebook is in the
``notebooks`` directory of the repository: run it from a clone of the repository with
sst1mpipe installed (``pip install -e .``) and the test data downloaded (``git lfs pull``).

Reading the raw data
--------------------

The SST-1M raw data (R0, ``SST1M<telescope>_<date>_<run>.fits.fz`` zfits files written by the
DigiCam camera server) are read with the ctapipe event source :class:`~sst1mpipe.io.sst1m_event_source.SST1MEventSource`
(see :doc:`../event_sources`). It is a ctapipe I/O plugin: `ctapipe.io.EventSource` opens the SST-1M
files by itself.

.. code-block:: python

   from ctapipe.io import EventSource

   with EventSource("SST1M2_20260120_1102.fits.fz", max_events=10) as source:
       print(source.target, source.pointing)
       for event in source:
           r0 = event.r0.tel[source.tel_id]
           print(event.index.event_id, event.trigger.event_type, r0.waveform.shape)

The notebook below goes through everything the event source gives: run information
(target, pointing, observation blocks, subarray, DigiCam configuration), the fields of the
events (index, trigger, R0 waveforms, baselines, trigger traces, monitoring), plots of the
waveforms and camera images, and the underlying ``protozfits`` data.

.. toctree::
   :maxdepth: 1

   sst1m_read_zfits_event_source

The DigiCam trigger
-------------------

The trigger traces recorded by the camera (trigger input traces of the patches, patch7 output), how the
patch7 trigger works, and how to emulate it from the raw data and from simulations with
``sst1mpipe.trigger`` (patch7, patch19, TDSCAN).

.. toctree::
   :maxdepth: 1

   sst1m_trigger_traces
