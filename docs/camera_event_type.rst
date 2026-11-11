.. _camera_event_type:

Trigger type of the camera events
=================================

The raw data (R0) of an event read by :class:`~sst1mpipe.io.sst1m_event_source.SST1MEventSource` (zfits files) or
:class:`~sst1mpipe.io.zmq_event_source.ZMQEventSource` (ZMQ stream) contain the event type sent by DigiCam, ``event.r0.tel[tel_id].event_type``. It is a
bit mask of the tpye trigger of the event, :class:`~sst1mpipe.io.containers.CameraEventType`.

.. autoclass:: sst1mpipe.io.containers.CameraEventType
   :no-members:
