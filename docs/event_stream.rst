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
- the messages published before the subscriber is connected are lost.

Messages
--------

One message per telescope event (topic ``tel<tel_id>``, e.g. ``tel21``) and one per stereo event
(topic ``stereo``). Each message has a header (msgpack) and the numpy arrays of the event without
copy. The data levels are ``trigger``, ``r0``, ``r1``, ``dl0``, ``dl1``, ``dl2`` and
``dl2_stereo``; their fields have the names of the columns of the ctapipe ``DataWriter``
(e.g. ``camera_frame_hillas_intensity``, ``<algorithm>_tel_energy``), in the units of the ctapipe
tables (deg, m, ns, TeV), the times in unix seconds.

To limit the bandwidth (the R0 waveforms of an event of DigiCam take about 130 kB), the full
events, with the waveforms and the images, are sent only every
``--ZMQEventPublisher.full_event_period`` seconds (1 s by default, 0 for all the events); the
other events have their trigger, DL1 parameters and DL2, enough for the statistics and the rates.

Receiving the events (subscriber)
--------------------

.. code-block:: python

    from sst1mpipe.io.zmq_event_stream import ZMQEventSubscriber

    with ZMQEventSubscriber("tcp://*:5561") as subscriber:
        while True:
            message = subscriber.receive(timeout=1.0)
            if message is not None:
                header, levels = message  # levels: {"dl1": {"hillas_intensity": ..., ...}, ...}

See ``sst1mpipe.io.zmq_event_stream`` for the details of the format.
