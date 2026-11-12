"""
Monitoring source reading data points from an OPC UA server (e.g. the drive of the telescope for the
real-time analysis), used by sst1mpipe-process
(``--monitoring-source OPCUAMonitoringSource --OPCUAMonitoringSource.host=<host> ...``).

The source is read only: it polls the values of the configured nodes at a fixed period in a
background thread, keeps them with the time at which they are read in a buffer, and fills the
monitoring containers of the events with the values read closest to the time of the event.

``asyncua`` is needed (``pip install sst1mpipe[opcua]``): it is imported when the source is created.
"""
import bisect
import threading
import time as clock
from collections import deque

import astropy.units as u
import numpy as np
from astropy.table import Table
from astropy.time import Time
from ctapipe.containers import TelescopeMonitoringContainer
from ctapipe.core.traits import AstroQuantity, Dict, Float, Int, Unicode
from ctapipe.io import MonitoringSource
from ctapipe.io.monitoringtypes import MonitoringType

__all__ = ["OPCUAMonitoringSource", "sercos_node_id"]

# monitoring type of the fields of the telescope monitoring container, by prefix of their path
MONITORING_TYPE_PREFIXES = {
    "pointing": MonitoringType.TELESCOPE_POINTINGS,
    "camera.pixel_statistics": MonitoringType.PIXEL_STATISTICS,
    "camera.coefficients": MonitoringType.CAMERA_COEFFICIENTS,
}


def sercos_node_id(bus, axis, parameter, namespace=7):
    """
    OPC UA node id of a Sercos parameter of a drive axis, e.g. ``sercos_node_id(0, 1, "S-0-0386")``
    (active position feedback of the axis 1 of the bus 0).
    """
    return f'ns={namespace};s=Sercos,{bus},{axis}.ParameterSet."{parameter}"'


def monitoring_type_of_field(field):
    """Monitoring type of a field of the telescope monitoring container (``pointing.altitude``), None if unknown"""
    for prefix, monitoring_type in MONITORING_TYPE_PREFIXES.items():
        if field == prefix or field.startswith(prefix + "."):
            return monitoring_type
    return None


def check_field(field):
    """Raise a ValueError if ``field`` is not a field of the telescope monitoring container"""
    container = TelescopeMonitoringContainer()
    *parents, name = field.split(".")
    for parent in parents:
        if parent not in container.fields:
            raise ValueError(f"{field!r}: {parent!r} is not a field of {type(container).__name__}")
        container = container[parent]
    if name not in container.fields:
        raise ValueError(f"{field!r}: {name!r} is not a field of {type(container).__name__}")


def set_field(container, field, value):
    """Set the field ``field`` (path ``a.b.c``) of a container"""
    *parents, name = field.split(".")
    for parent in parents:
        container = container[parent]
    container[name] = value


class OPCUAMonitoringSource(MonitoringSource):
    """
    Monitoring data read from an OPC UA server (read only), e.g. the position of the drive of the
    telescope in the real-time analysis.

    The values of the nodes ``nodes`` (field of the telescope monitoring container
    ``event.monitoring.tel[tel_id]`` -> OPC UA node id) are read every ``polling_period`` seconds
    in a background thread, all in one request, and kept with the time at which they are read
    (time of the computer) for ``buffer_duration`` seconds. The events are filled with the values
    read closest to their trigger time, within ``timestamp_tolerance``; otherwise (e.g. server not
    reachable, or a node not readable) the field is not filled (default value, NaN for the pointing).
    The connection is retried every ``reconnect_period`` seconds if it is lost.

    A value is converted as ``value * scales[field] + offsets[field]``, in ``units[field]`` if given
    (e.g. the positions of the drive in deg). Example of configuration (json)::

        "OPCUAMonitoringSource": {
            "host": "drive.sst1m.local",
            "port": 24840,
            "polling_period": 1.0,
            "nodes": {
                "pointing.altitude": "ns=7;s=Sercos,0,1.ParameterSet.\\"S-0-0386\\"",
                "pointing.azimuth": "ns=7;s=Sercos,0,2.ParameterSet.\\"S-0-0386\\""
            },
            "units": {"pointing.altitude": "deg", "pointing.azimuth": "deg"}
        }

    The events are not skipped by sst1mpipe-process when the values are not available
    (``has_pointings`` is False): the real-time analysis goes on without the server.
    """

    host = Unicode("localhost", help="Host name or address of the OPC UA server").tag(config=True)
    port = Int(4840, help="Port of the OPC UA server").tag(config=True)
    nodes = Dict(
        key_trait=Unicode(),
        value_trait=Unicode(),
        default_value={},
        help=(
            "OPC UA node ids of the fields of the telescope monitoring container, by path of the field"
            " (e.g. {'pointing.altitude': 'ns=7;s=Sercos,0,1.ParameterSet.\"S-0-0386\"'})"
        ),
    ).tag(config=True)
    units = Dict(
        key_trait=Unicode(), value_trait=Unicode(), default_value={},
        help="Units of the values of the fields (e.g. {'pointing.altitude': 'deg'}), without unit if not given",
    ).tag(config=True)
    scales = Dict(
        key_trait=Unicode(), value_trait=Float(), default_value={},
        help="Scale factors of the values of the fields (value * scale + offset), 1 if not given",
    ).tag(config=True)
    offsets = Dict(
        key_trait=Unicode(), value_trait=Float(), default_value={},
        help="Offsets of the values of the fields (value * scale + offset), 0 if not given",
    ).tag(config=True)
    tel_id = Int(
        default_value=None, allow_none=True,
        help="Telescope whose monitoring container is filled, the telescope of the subarray if None",
    ).tag(config=True)
    polling_period = Float(1.0, min=0.01, help="Period of the reading of the values (s)").tag(config=True)
    timeout = Float(4.0, min=0.1, help="Timeout of the connection and of the requests to the server (s)").tag(config=True)
    reconnect_period = Float(10.0, min=0.1, help="Period of the attempts to connect to the server (s)").tag(config=True)
    buffer_duration = Float(
        600.0, min=0,
        help="Duration for which the values read are kept (s), to fill events processed later than their time",
    ).tag(config=True)
    timestamp_tolerance = AstroQuantity(
        default_value=2.0 * u.s,
        physical_type=u.physical.time,
        help="Maximum time difference between an event and the closest values read for them to be used",
    ).tag(config=True)

    def __init__(self, subarray=None, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        try:
            from asyncua.sync import Client  # noqa: F401
        except ImportError as error:
            raise ImportError("OPCUAMonitoringSource needs asyncua: pip install sst1mpipe[opcua]") from error
        for field in self.nodes:
            check_field(field)
        for name in ("units", "scales", "offsets"):
            unknown = set(getattr(self, name)) - set(self.nodes)
            if unknown:
                raise ValueError(f"OPCUAMonitoringSource.{name}: fields without node: {sorted(unknown)}")
        self._units = {field: u.Unit(unit) for field, unit in self.units.items()}
        if self.tel_id is None:
            if self.subarray is None or len(self.subarray.tel_ids) != 1:
                raise ValueError("OPCUAMonitoringSource.tel_id must be given for a subarray of several telescopes")
            self._tel_id = int(self.subarray.tel_ids[0])
        else:
            self._tel_id = self.tel_id
        if self.timestamp_tolerance.to_value(u.s) < self.polling_period:
            self.log.warning(
                "OPCUAMonitoringSource.timestamp_tolerance (%s) is shorter than the polling period (%s s):"
                " events may not be filled", self.timestamp_tolerance, self.polling_period,
            )
        self._fields = list(self.nodes)
        # samples (unix time of the reading, values of the fields), in time order
        self._times = deque()
        self._samples = deque()
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._first_poll = threading.Event()
        self._client = None
        self._unreadable = set()
        self._thread = None
        if self._fields:
            self._thread = threading.Thread(target=self._poll_loop, name="OPCUAMonitoringSource", daemon=True)
            self._thread.start()
            # values of the first events, if the server answers
            self._first_poll.wait(self.timeout + self.polling_period)
        else:
            self.log.warning("OPCUAMonitoringSource without nodes: no monitoring data")

    @property
    def url(self):
        """URL of the OPC UA server"""
        return f"opc.tcp://{self.host}:{self.port}"

    @property
    def monitoring_types(self):
        types = {monitoring_type_of_field(field) for field in self._fields}
        return tuple(t for t in MonitoringType if t in types)

    @property
    def is_connected(self):
        """True if the source is connected to the server"""
        return self._client is not None

    # ---------------------------------------------------------------- reading of the server

    def _connect(self):
        from asyncua.sync import Client

        # the event loop of the source is kept between the connections
        client = Client(self.url, timeout=self.timeout, tloop=self._tloop)
        try:
            client.connect()
        except Exception:
            client.disconnect_socket()
            raise
        try:
            nodes = [client.get_node(self.nodes[field]) for field in self._fields]
        except Exception:
            self._client = client
            self._disconnect()
            raise
        self._client = client
        self._opcua_nodes = nodes
        self.log.info("Connected to the OPC UA server %s", self.url)

    def _disconnect(self):
        client, self._client = self._client, None
        if client is not None:
            try:
                client.disconnect()
            except Exception as error:  # the connection may already be lost
                self.log.debug("Disconnection from %s: %s", self.url, error)

    def _read(self):
        """Values of the fields (NaN if not readable), and unix time of the reading"""
        start = clock.time()
        data_values = self._client.read_attributes(self._opcua_nodes)
        read_time = 0.5 * (start + clock.time())
        values = []
        for field, data_value in zip(self._fields, data_values, strict=True):
            if data_value.StatusCode is not None and not data_value.StatusCode.is_good():
                if field not in self._unreadable:
                    self.log.warning(
                        "OPC UA node %s (%s) not readable: %s", self.nodes[field], field, data_value.StatusCode.name,
                    )
                    self._unreadable.add(field)
                values.append(np.nan)
                continue
            self._unreadable.discard(field)
            values.append(self._convert(field, data_value.Value.Value))
        return read_time, values

    def _convert(self, field, value):
        if field in self.scales or field in self.offsets or field in self._units:
            value = float(value) * self.scales.get(field, 1.0) + self.offsets.get(field, 0.0)
        if field in self._units:
            value = u.Quantity(value, self._units[field])
        return value

    def _poll_loop(self):
        from asyncua.sync import ThreadLoop

        self._tloop = ThreadLoop(timeout=self.timeout + 1)
        self._tloop.start()
        try:
            self._poll()
        finally:
            self._disconnect()
            self._tloop.stop()
            self._first_poll.set()

    def _poll(self):
        next_connection = 0.0
        connection_failed = False
        while not self._stop.is_set():
            start = clock.monotonic()
            if self._client is None and start >= next_connection:
                try:
                    self._connect()
                    connection_failed = False
                except Exception as error:
                    # warning at the first failure only, while the server is not reachable
                    log = self.log.debug if connection_failed else self.log.warning
                    log(
                        "Cannot connect to the OPC UA server %s (%s: %s), retry every %s s",
                        self.url, type(error).__name__, error, self.reconnect_period,
                    )
                    connection_failed = True
                    next_connection = start + self.reconnect_period
            if self._client is not None:
                try:
                    self._add_sample(*self._read())
                except Exception as error:
                    self.log.warning(
                        "Reading of the OPC UA server %s failed (%s: %s), reconnection",
                        self.url, type(error).__name__, error,
                    )
                    self._disconnect()
            self._first_poll.set()
            self._stop.wait(max(0.0, self.polling_period - (clock.monotonic() - start)))

    def _add_sample(self, read_time, values):
        with self._lock:
            self._times.append(read_time)
            self._samples.append(values)
            while self._times and self._times[0] < read_time - self.buffer_duration:
                self._times.popleft()
                self._samples.popleft()

    # ---------------------------------------------------------------- access to the values

    def _check_tel_id(self, tel_id):
        if tel_id is not None and tel_id != self._tel_id:
            raise KeyError(f"OPCUAMonitoringSource: no monitoring data for tel_id {tel_id} (tel_id {self._tel_id})")

    def _type_fields(self, monitoring_type):
        if monitoring_type not in self.monitoring_types:
            raise KeyError(f"OPCUAMonitoringSource: no monitoring data of type {monitoring_type}")
        return [i for i, field in enumerate(self._fields) if monitoring_type_of_field(field) == monitoring_type]

    def values_at(self, time):
        """
        Values of all the fields read closest to ``time`` (scalar `~astropy.time.Time`), by field;
        empty if no values were read within ``timestamp_tolerance``
        """
        unix = time.unix
        with self._lock:
            index = bisect.bisect_left(self._times, unix)
            candidates = [i for i in (index - 1, index) if 0 <= i < len(self._times)]
            if not candidates:
                return {}
            closest = min(candidates, key=lambda i: abs(self._times[i] - unix))
            if abs(self._times[closest] - unix) > self.timestamp_tolerance.to_value(u.s):
                return {}
            return dict(zip(self._fields, self._samples[closest], strict=True))

    def get_table(self, monitoring_type, tel_id=None, **kwargs):
        """Values read and kept in the buffer of the fields of a monitoring type, with their time"""
        self._check_tel_id(tel_id)
        indices = self._type_fields(monitoring_type)
        with self._lock:
            times, samples = list(self._times), list(self._samples)
        table = Table()
        table["time"] = Time(times, format="unix")
        for i in indices:
            column = [sample[i] for sample in samples]
            unit = self._units.get(self._fields[i])
            if unit is not None:
                column = u.Quantity([c.to_value(unit) if isinstance(c, u.Quantity) else c for c in column], unit)
            table[self._fields[i]] = column
        return table

    def get_values(self, monitoring_type, time, tel_id=None, **kwargs):
        """
        Values of the fields of a monitoring type (by name of the field in its container, e.g.
        ``altitude``) read closest to ``time``. Missing values (no values within ``timestamp_tolerance``) are NaN.
        """
        self._check_tel_id(tel_id)
        indices = self._type_fields(monitoring_type)
        times = time.reshape(-1)
        rows = [self.values_at(t) for t in times]
        result = {}
        for i in indices:
            field = self._fields[i]
            name = field.split(".")[-1]
            unit = self._units.get(field)
            values = [row.get(field, np.nan) for row in rows]
            if unit is not None:
                values = u.Quantity([v if isinstance(v, u.Quantity) else v * unit for v in values], unit)
            result[name] = values[0] if time.isscalar else values
        return result

    def fill_monitoring_container(self, event):
        """Fill the monitoring container of the telescope with the values read closest to the time of the event"""
        values = self.values_at(event.trigger.time)
        monitoring = event.monitoring.tel[self._tel_id]
        for field, value in values.items():
            if isinstance(value, float) and np.isnan(value):
                continue
            set_field(monitoring, field, value)

    def close(self):
        """Stop the reading of the server and disconnect"""
        self._stop.set()
        if self._thread is not None and self._thread.is_alive():
            self._thread.join(self.timeout + self.polling_period + 1)
