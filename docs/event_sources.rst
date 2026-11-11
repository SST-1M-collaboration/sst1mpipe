Event sources
=============

The SST-1M data are read with two ctapipe event sources of sst1mpipe, registered as ctapipe I/O
plugins: `ctapipe.io.EventSource` chooses the right one by itself from the input.

* :class:`~sst1mpipe.io.sst1m_event_source.SST1MEventSource` reads the raw data files (R0,
  ``SST1M<telescope>_<date>_<run>.fits.fz``).
* :class:`~sst1mpipe.io.zmq_event_source.ZMQEventSource` reads the events sent in a ZMQ stream, e.g.
  by the camera server for the real time analysis: DigiCam camera events (R0, as in the files) and
  CTAO R1 and DL0 events.

Both fill the R0 data of the DigiCam camera events in the same way: see the tutorial
:doc:`tutorials/sst1m_read_zfits_event_source` for what they contain.

.. code-block:: python

   from ctapipe.io import EventSource

   with EventSource("SST1M2_20260120_1102.fits.fz", max_events=10) as source:
       for event in source:
           print(event.index.event_id, event.trigger.event_type)

SST1MEventSource
----------------

.. autoclass:: sst1mpipe.io.sst1m_event_source.SST1MEventSource
   :members:
   :exclude-members: get_array_event

ZMQEventSource
--------------

.. autoclass:: sst1mpipe.io.zmq_event_source.ZMQEventSource
   :members:
