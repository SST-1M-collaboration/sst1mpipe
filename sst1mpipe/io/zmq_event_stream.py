"""
Stream of the processed events of the pipeline (ZMQ), e.g. for the viewer of the events.

The pipeline publishes its events (``ZMQEventPublisher``, PUB socket which *connects*) and a single
subscriber gets them (``ZMQEventSubscriber``, SUB socket which *binds*): several processes of the
pipeline (e.g. the workers of a telescope, which share its camera stream in PULL) can publish to
the same subscriber, each event once. A publisher never waits for the subscriber: when the subscriber
is absent or too slow, the messages are dropped (the subscriber counts them, see below).

Messages
--------
One message per telescope event and one per stereo event, multipart:

1. the topic: ``tel<tel_id>`` (e.g. ``tel21``) or ``stereo``;
2. the header (msgpack): format ``version``, ``source`` (name of the publisher), ``count``
   (number of the message of the publisher for this topic: a gap is a lost message), ``obs_id``, ``event_id``,
   ``tel_id`` (0 for stereo), and for each data level of the event its ``scalars``
   ``{name: value}`` and the descriptions of its ``arrays`` ``[name, dtype, shape, frame]``
   (``frame``: index of the frame of the array among the array frames: an array of several
   levels, e.g. the DL0 waveform equal to the R1 waveform, is sent once);
3. one frame per array, the raw bytes of the numpy array (no copy).

The data levels are ``trigger``, ``r0``, ``r1``, ``dl0``, ``dl1``, ``dl2`` (of the telescope)
and ``dl2_stereo``. Their fields have the names of the columns of the ctapipe ``DataWriter``
(e.g. ``hillas_intensity``, ``camera_frame_hillas_intensity``, ``<algorithm>_tel_energy``,
``<algorithm>_alt`` for stereo), with the quantities in the units of the ctapipe tables (deg, m,
ns, TeV) and the times in unix seconds.

To limit the bandwidth, the publisher sends the full events (waveforms, images) only every
``full_event_period`` seconds: the first event after the period with a DL1 image (a finite Hillas
intensity), or any event after two periods. The other events have their trigger, DL1 parameters and
DL2 only (enough for the statistics and the rates).
"""

import os
import socket
import time
from enum import Enum

import astropy.units as u
import msgpack
import numpy as np
import zmq
from astropy.time import Time
from ctapipe.core import Component, Container
from ctapipe.core.traits import Float, Int, Unicode

__all__ = [
    "FORMAT_VERSION", "STEREO_TOPIC", "LEVELS", "ZMQEventPublisher", "ZMQEventSubscriber", "event_messages",
    "encode", "decode", "telescope_topic",
]

#: version of the format of the messages
FORMAT_VERSION = 1
#: topic of the stereo events
STEREO_TOPIC = "stereo"
#: data levels of the messages
LEVELS = ("trigger", "r0", "r1", "dl0", "dl1", "dl2", "dl2_stereo")
#: levels sent only in the full events
FULL_LEVELS = ("r0", "r1", "dl0")
#: DL1 fields sent only in the full events (the images)
IMAGE_FIELDS = ("image", "peak_time", "image_mask")

#: units of the values of the quantities, by physical type (those of the ctapipe tables)
UNITS = {"angle": u.deg, "length": u.m, "time": u.ns, "energy": u.TeV}
_UNIT_CACHE = {}


def telescope_topic(tel_id):
    """Topic of the messages of a telescope, e.g. ``tel21``"""
    return f"tel{int(tel_id)}"


def _unit_factor(unit):
    """
    Factor converting the values of a quantity of ``unit`` to the unit of its physical type
    (`UNITS`), None to keep them (cached: the physical type and the conversion are slow to get)
    """
    factor = _UNIT_CACHE.get(unit, False)
    if factor is False:
        target = UNITS.get(str(unit.physical_type))
        factor = _UNIT_CACHE[unit] = None if target is None or unit == target else unit.to(target)
    return factor


#: julian date of the unix epoch
_UNIX_EPOCH_JD = 2440587.5
#: day (unix_tai // 86400) -> unix - unix_tai (the leap seconds change only at the start of a day)
_LEAP_SECONDS = {}
#: plain types sent as they are
_PLAIN_TYPES = (bool, int, float, str)


def _unix(value):
    """
    Unix time (s) of a scalar `~astropy.time.Time`: its TAI unix time (fast) corrected by the leap
    seconds of its day (cached; the UTC conversion of astropy, ~70 us, is done once per day)
    """
    if value.scale != "tai" or not value.isscalar:
        return value.unix
    # the julian date of the internal representation: 1 us (the public jd1 / jd2 / unix_tai: 7 to 16 us)
    jd = getattr(value, "_time", value)
    unix_tai = float((jd.jd1 - _UNIX_EPOCH_JD + jd.jd2) * 86400.0)
    day = unix_tai // 86400.0
    offset = _LEAP_SECONDS.get(day)
    if offset is None:
        offset = _LEAP_SECONDS[day] = float(value.unix) - unix_tai
    return unix_tai + offset


def _value(value):
    """Plain value (number, string, numpy array) of a field value, None if it cannot be sent"""
    kind = type(value)
    if kind in _PLAIN_TYPES:  # the most frequent (fast path)
        return value
    if isinstance(value, np.generic):  # numpy scalar
        return value.item() if value.dtype != object else None
    if value is None or isinstance(value, (Container, dict)):
        return None
    if isinstance(value, u.Quantity):
        factor = _unit_factor(value.unit)
        value = value.value if factor is None else value.value * factor
    elif isinstance(value, Time):
        return _unix(value)  # times in unix seconds
    elif isinstance(value, Enum):
        value = value.value
    if isinstance(value, (bool, int, float, str, np.bool_, np.integer, np.floating)):
        return value.item() if isinstance(value, np.generic) else value
    array = np.asarray(value)
    if array.dtype == object:
        return None
    return array.item() if array.ndim == 0 else array


def _fields(container, recursive=False):
    """Fields ``{name: value}`` of a container (flattened with the prefixes, as by the DataWriter)"""
    items = container.as_dict(recursive=True, flatten=True, add_prefix=True).items() if recursive \
        else container.items()
    return {name: value for name, value in ((n, _value(v)) for n, v in items) if value is not None}


def _dl2_row(reconstructed, stereo):
    """Fields of the DL2 of an event: ``<algorithm>_tel_<field>`` per telescope, ``<algorithm>_<field>`` stereo"""
    row = {}
    for group, algorithm_map in reconstructed.items():
        # the fields of the impact have a prefix, as in the ctapipe tables
        infix = "impact_" if group == "impact" else ""
        for algorithm, container in algorithm_map.items():
            prefix = f"{algorithm}_" if stereo else f"{algorithm}_tel_"
            for field, value in container.items():
                value = _value(value)
                if value is not None:
                    row[f"{prefix}{infix}{field}"] = value
    return row


def event_messages(event, full=True, tel_ids=None, stereo=True):
    """
    The contents of the messages of an event.

    Parameters
    ----------
    event : ctapipe.containers.ArrayEventContainer
        The processed event (e.g. a `~sst1mpipe.io.containers.SST1MArrayEventContainer`).
    full : bool
        Send the waveforms (R0, R1, DL0) and the DL1 images; else only the trigger, the DL1
        parameters and the DL2.
    tel_ids : iterable or None
        The telescopes sent, all if None.
    stereo : bool
        Send the stereo DL2.

    Returns
    -------
    list of (topic, tel_id, levels)
        ``levels`` is ``{level: {field: value}}``.
    """
    tel_ids = None if tel_ids is None else {int(t) for t in tel_ids}
    telescopes = sorted({*event.trigger.tel, *event.r0.tel, *event.r1.tel, *event.dl0.tel, *event.dl1.tel,
                         *event.dl2.tel})
    event_type = _value(event.trigger.event_type)
    messages = []
    for tel_id in telescopes:
        if tel_ids is not None and tel_id not in tel_ids:
            continue
        levels = {}
        if tel_id in event.trigger.tel:
            trigger = _fields(event.trigger.tel[tel_id])
            trigger["event_type"] = event_type
            levels["trigger"] = trigger
        if full:
            for level in FULL_LEVELS:
                container = getattr(event, level).tel.get(tel_id)
                if container is not None:
                    fields = _fields(container)
                    if fields:
                        levels[level] = fields
            _share_waveform(levels)
        if tel_id in event.dl1.tel:
            dl1 = event.dl1.tel[tel_id]
            fields = _fields(dl1.parameters, recursive=True) if dl1.parameters is not None else {}
            for name in ("is_valid", *(IMAGE_FIELDS if full else ())):
                value = _value(getattr(dl1, name)) if name in dl1.fields else None
                if value is not None:
                    fields[name] = value
            levels["dl1"] = fields
        if tel_id in event.dl2.tel:
            row = _dl2_row(event.dl2.tel[tel_id], stereo=False)
            if row:
                levels["dl2"] = row
        if levels:
            messages.append((telescope_topic(tel_id), tel_id, levels))
    if stereo:
        row = _dl2_row(event.dl2.stereo, stereo=True)
        if row:
            row["time"] = _value(event.trigger.time)
            messages.append((STEREO_TOPIC, 0, {"dl2_stereo": row}))
    return messages


def _share_waveform(levels):
    """
    The DL0 waveform equal to the R1 waveform (no data volume reduction) is the same array, sent
    once by `encode`
    """
    r1, dl0 = levels.get("r1", {}).get("waveform"), levels.get("dl0", {}).get("waveform")
    if isinstance(r1, np.ndarray) and isinstance(dl0, np.ndarray) and r1 is not dl0 and np.array_equal(r1, dl0):
        levels["dl0"]["waveform"] = r1


def _has_image(event):
    """Whether a telescope of the event has a DL1 image (a finite Hillas intensity)"""
    for dl1 in event.dl1.tel.values():
        parameters = dl1.parameters
        intensity = parameters.hillas.intensity if parameters is not None else None
        if intensity is not None and np.isfinite(intensity):
            return True
    return False


def encode(topic, levels, obs_id, event_id, tel_id, source="", count=0):
    """
    Frames of a message (see the module).

    Parameters
    ----------
    topic : str
    levels : dict
        ``{level: {field: value}}``, the values being numbers, strings or numpy arrays.
    obs_id, event_id, tel_id : int
    source : str
        Name of the publisher.
    count : int
        Number of the message of the publisher for this topic.

    Returns
    -------
    list
        The frames: topic, header, then the arrays (sent without copy).
    """
    arrays, frames = [], {}  # frames: id of an array -> index of its frame
    header = dict(version=FORMAT_VERSION, source=source, count=int(count), obs_id=int(obs_id),
                  event_id=int(event_id), tel_id=int(tel_id), levels={})
    for level, fields in levels.items():
        scalars, descriptions = {}, []
        for name, value in fields.items():
            if isinstance(value, np.ndarray):
                index = frames.get(id(value))
                if index is None:  # an array of several levels is sent once
                    index = frames[id(value)] = len(arrays)
                    arrays.append(np.ascontiguousarray(value))
                descriptions.append([name, value.dtype.str, list(value.shape), index])
            else:
                scalars[name] = value
        header["levels"][level] = dict(scalars=scalars, arrays=descriptions)
    return [topic.encode(), msgpack.packb(header, use_bin_type=True), *arrays]


def decode(frames):
    """
    Header and data levels of a message (`encode`).

    Returns
    -------
    header : dict
        The header (``source``, ``count``, ``obs_id``, ``event_id``, ``tel_id``, ...).
    levels : dict
        ``{level: row}``, each row having its fields, ``obs_id``, ``event_id`` and ``tel_id``.
    """
    header = msgpack.unpackb(bytes(frames[1]), raw=False)
    if header.get("version") != FORMAT_VERSION:
        raise ValueError(f"Unknown format version {header.get('version')} (known: {FORMAT_VERSION})")
    keys = dict(obs_id=header["obs_id"], event_id=header["event_id"], tel_id=header["tel_id"])
    buffers = frames[2:]
    arrays = {}  # index of the frame -> array (an array of several levels is the same array)
    levels = {}
    for level, content in header["levels"].items():
        row = dict(content["scalars"])
        for name, dtype, shape, index in content["arrays"]:
            array = arrays.get(index)
            if array is None:
                frame = buffers[index]
                data = frame.buffer if isinstance(frame, zmq.Frame) else frame
                array = arrays[index] = np.frombuffer(data, dtype=np.dtype(dtype)).reshape(shape)
            row[name] = array
        row.update(keys)
        levels[level] = row
    return header, levels


class ZMQEventPublisher(Component):
    """
    Publish the processed events in a ZMQ stream (PUB socket connecting to ``url``): see the
    module. Call it with each event; without ``url`` nothing is sent.
    """

    url = Unicode(
        default_value=None, allow_none=True,
        help="ZMQ endpoint of the subscriber of the events (e.g. tcp://viewer-host:5561); none: no stream",
    ).tag(config=True)
    full_event_period = Float(
        default_value=1.0,
        help="period (s) of the full events, with the waveforms and the images (0: all the events)",
    ).tag(config=True)
    send_high_water_mark = Int(
        default_value=1000, help="maximum number of messages queued for the subscriber (the others are dropped)",
    ).tag(config=True)
    source_name = Unicode(
        default_value=f"{socket.gethostname()}:{os.getpid()}", help="name of this publisher in the messages",
    ).tag(config=True)

    def __init__(self, config=None, parent=None, **kwargs):
        super().__init__(config=config, parent=parent, **kwargs)
        self.counts = {}  # topic -> number of messages sent
        self._last_full = -np.inf
        self._context = self._socket = None
        if self.url:
            self._context = zmq.Context.instance()
            self._socket = self._context.socket(zmq.PUB)
            self._socket.setsockopt(zmq.SNDHWM, self.send_high_water_mark)
            self._socket.setsockopt(zmq.LINGER, 0)
            self._socket.connect(self.url)
            self.log.info("Publishing the events to %s", self.url)

    @property
    def enabled(self):
        """Whether the events are published"""
        return self._socket is not None

    def __call__(self, event, tel_ids=None, stereo=True):
        """Publish an event (its telescopes ``tel_ids``, all if None, and its stereo DL2)"""
        if self._socket is None:
            return
        now = time.monotonic()
        # a full event after the period: one with an image, or any after two periods
        elapsed = now - self._last_full
        full = elapsed >= self.full_event_period and (elapsed >= 2 * self.full_event_period or _has_image(event))
        if full:
            self._last_full = now
        obs_id, event_id = event.index.obs_id, event.index.event_id
        for topic, tel_id, levels in event_messages(event, full=full, tel_ids=tel_ids, stereo=stereo):
            count = self.counts.get(topic, 0)
            frames = encode(topic, levels, obs_id, event_id, tel_id, self.source_name, count)
            self.counts[topic] = count + 1
            # PUB never blocks: the messages which do not fit in the queue are dropped
            self._socket.send_multipart(frames, copy=False)

    def close(self):
        """Close the socket"""
        if self._socket is not None:
            self._socket.close(linger=0)
            self._socket = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class ZMQEventSubscriber:
    """
    Receive the events published by `ZMQEventPublisher` (SUB socket binding ``url``, e.g.
    ``tcp://*:5561``): several publishers can connect to it.

    Parameters
    ----------
    url : str
        Endpoint bound by the subscriber.
    topics : iterable of str
        Topics received (e.g. ``["tel21"]``), all if empty.
    receive_high_water_mark : int
        Maximum number of messages queued (the others are dropped).

    Attributes
    ----------
    received : int
        Number of messages received.
    lost : dict
        ``{source: number of messages lost}``, from the gaps of the message counts of each
        publisher and topic (messages dropped because the subscriber was too slow).
    """

    def __init__(self, url, topics=(), receive_high_water_mark=10000):
        self.url = url
        self._context = zmq.Context.instance()
        self.socket = self._context.socket(zmq.SUB)
        self.socket.setsockopt(zmq.RCVHWM, receive_high_water_mark)
        self.socket.setsockopt(zmq.LINGER, 0)
        for topic in topics or ("",):
            self.socket.setsockopt(zmq.SUBSCRIBE, topic.encode())
        self.socket.bind(url)
        #: the endpoint bound (e.g. with the port chosen by the system for ``tcp://127.0.0.1:*``)
        self.endpoint = self.socket.getsockopt_string(zmq.LAST_ENDPOINT)
        self.received = 0
        self.lost = {}
        self._next_count = {}

    def receive(self, timeout=None):
        """
        Receive a message.

        Parameters
        ----------
        timeout : float or None
            Maximum waiting time (s), None to wait without limit.

        Returns
        -------
        (header, levels) or None
            The decoded message (`decode`), None if nothing was received in time.
        """
        try:  # a message waiting: no poll (one system call less per message)
            frames = self.socket.recv_multipart(flags=zmq.NOBLOCK, copy=False)
        except zmq.Again:
            if timeout is not None and not self.socket.poll(int(timeout * 1000)):
                return None
            frames = self.socket.recv_multipart(copy=False)
        header, levels = decode(frames)
        self.received += 1
        source, count = header["source"], header["count"]
        key = (source, bytes(frames[0]))  # the messages are numbered per publisher and topic
        expected = self._next_count.get(key)
        if expected is not None and count > expected:
            self.lost[source] = self.lost.get(source, 0) + count - expected
        self._next_count[key] = count + 1
        return header, levels

    def close(self):
        """Close the socket"""
        self.socket.close(linger=0)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
