Stream of the processed events
==============================

``sst1mpipe-process`` can publish the events it processes in a ZMQ stream, e.g. for the viewer
of the events (``sst1mpipe-viewer``), with ``--publish``:

.. code-block:: console

    sst1mpipe-process -i tcp://camera-server:24593 -o events.dl1.h5 --config sst1mpipe_rta_config.json \
        --publish tcp://viewer-host:5561

The subscriber (e.g. the viewer) *binds* the endpoint and the processes of the pipeline *connect*
to it (PUB / SUB sockets of ZMQ):

- several processes can publish to the same subscriber, e.g. the workers of a telescope sharing its
  camera stream (PULL): the subscriber gets each event once, from the worker which processed it.
  The events of different workers can arrive in a different order than their timestamps;
- a publisher never waits for the subscriber: if the subscriber is absent or too slow, the messages
  are dropped. The messages of each publisher are numbered, so that the subscriber counts the
  messages lost;
- the messages published before the subscriber is connected are lost;
- when a publisher is closed, the messages still queued are sent during at most
  ``--ZMQEventPublisher.close_linger`` seconds (1 s by default).

The publisher is ``sst1mpipe.io.zmq_event_stream.ZMQEventPublisher`` (a ctapipe component, called
with each event), the subscriber ``sst1mpipe.io.zmq_event_stream.ZMQEventSubscriber``.

Messages
--------

One message per telescope event (topic ``tel<tel_id>``, e.g. ``tel21``) and one per stereo event
(topic ``stereo``). A message is multipart:

1. the topic;
2. the header (msgpack): format ``version``, ``source`` (name of the publisher, ``<host>:<pid>`` by
   default), ``count`` (number of the message of the publisher for this topic: a gap is a lost
   message), ``obs_id``, ``event_id``, ``tel_id`` (0 for stereo), and for each data level of the
   event its ``scalars`` ``{name: value}`` and the descriptions of its ``arrays``
   ``[name, dtype, shape, frame]``. ``frame`` is the index of the frame of the array among the
   array frames: an array of several levels, e.g. the DL0 waveform equal to the R1 waveform, is
   sent once;
3. one frame per array: the raw bytes of the numpy array, sent without copy.

The data levels are ``trigger``, ``r0``, ``r1``, ``dl0``, ``dl1``, ``dl2`` (of the telescope) and
``dl2_stereo``. Their fields have the names of the columns of the ctapipe ``DataWriter``
(e.g. ``hillas_intensity``, ``camera_frame_hillas_intensity``, ``<algorithm>_tel_energy``,
``<algorithm>_alt`` for stereo), with the quantities in the units of the ctapipe tables (deg, m,
ns, TeV) and the times in unix seconds. The fields which are ``None`` (e.g. trigger traces not
sent by the camera) are not in the messages.

Full events
-----------

To limit the bandwidth (a full event of DigiCam takes about 770 kB: R0, R1 and DL0 waveforms,
trigger traces, DL1 images), the full events are sent only every
``--ZMQEventPublisher.full_event_period`` seconds (1 s by default, 0 for all the events): the first
event after the period with a DL1 image (a finite Hillas intensity), or any event after two
periods. The other events have their trigger, DL1 parameters and DL2 (about 1.5 kB), enough for
the statistics and the rates.

Receiving the events (subscriber)
---------------------------------

.. code-block:: python

    from sst1mpipe.io.zmq_event_stream import ZMQEventSubscriber

    with ZMQEventSubscriber("tcp://*:5561") as subscriber:
        while True:
            message = subscriber.receive(timeout=1.0)
            if message is not None:
                header, levels = message  # levels: {"dl1": {"hillas_intensity": ..., ...}, ...}

