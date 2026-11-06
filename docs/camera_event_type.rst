.. _camera_event_type:

Trigger type of the camera events
=================================

The raw data (R0) of an event read by ``SST1MEventSource`` (zfits files) or ``ZMQEventSource``
(ZMQ stream) contain the event type sent by DigiCam, ``event.r0.tel[tel_id].event_type``. It is a
bit mask of the tpye trigger of the event, :class:`~sst1mpipe.io.containers.CameraEventType`.


.. code-block:: python

    from sst1mpipe.io.containers import CameraEventType
    from sst1mpipe.io.sst1m_event_source import SST1MEventSource

    with SST1MEventSource("SST1M1_20260121_0206.fits.fz") as source:
        for event in source:
            event_type = event.r0.tel[21].event_type
            if CameraEventType.MUON_TRIGGER in event_type:
                ...

.. autoclass:: sst1mpipe.io.containers.CameraEventType
   :no-members:
