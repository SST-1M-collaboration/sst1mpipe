"""
OPCUAMonitoringSource: monitoring data read from a local OPC UA server mimicking the drive
(Sercos parameters of the axes)
"""
import socket
import time

import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time
from ctapipe.containers import ArrayEventContainer
from ctapipe.io import MonitoringSource
from ctapipe.io.monitoringtypes import MonitoringType

from sst1mpipe.io.opcua_monitoring_source import OPCUAMonitoringSource, sercos_node_id
from sst1mpipe.io.sst1m_event_source import SUBARRAY_DESCRIPTION

asyncua_sync = pytest.importorskip("asyncua.sync")
ua = pytest.importorskip("asyncua.ua")

TEL_ID = 21
SUBARRAY = SUBARRAY_DESCRIPTION.select_subarray([TEL_ID])
NAMESPACE_URI = "urn:sst1mpipe:test:drive"
POLLING_PERIOD = 0.1


def free_port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


class DriveServer:
    """OPC UA server with the active position feedback (S-0-0386) of the two axes"""

    def __init__(self, port):
        self.port = port
        self.server = asyncua_sync.Server()
        self.server.set_endpoint(f"opc.tcp://127.0.0.1:{port}")
        self.namespace = self.server.register_namespace(NAMESPACE_URI)
        objects = self.server.nodes.objects
        self.variables = {}
        for axis, value in ((1, 70.0), (2, 180.0)):
            node_id = ua.NodeId(f'Sercos,0,{axis}.ParameterSet."S-0-0386"', self.namespace)
            self.variables[axis] = objects.add_variable(node_id, f"axis_{axis}_position", value)
        # node with a text value
        objects.add_variable(ua.NodeId("status", self.namespace), "status", "tracking")
        # the port of a stopped server may take a moment to be free again
        for attempt in range(20):
            try:
                self.server.start()
                break
            except OSError:
                if attempt == 19:
                    raise
                time.sleep(0.25)

    def node_id(self, axis, parameter="S-0-0386"):
        return sercos_node_id(0, axis, parameter, namespace=self.namespace)

    def set(self, axis, value):
        self.variables[axis].write_value(value)

    def stop(self):
        self.server.stop()


@pytest.fixture(scope="module")
def server():
    drive = DriveServer(free_port())
    yield drive
    drive.stop()


def pointing_config(server, **kwargs):
    return dict(
        port=server.port,
        host="127.0.0.1",
        polling_period=POLLING_PERIOD,
        timestamp_tolerance=1 * u.s,
        nodes={"pointing.altitude": server.node_id(1), "pointing.azimuth": server.node_id(2)},
        units={"pointing.altitude": "deg", "pointing.azimuth": "deg"},
        **kwargs,
    )


def event_at(time):
    event = ArrayEventContainer()
    event.trigger.time = time
    return event


def wait_for(condition, timeout=5.0):
    start = time.monotonic()
    while time.monotonic() - start < timeout:
        if condition():
            return True
        time.sleep(0.05)
    return False


def test_sercos_node_id():
    assert sercos_node_id(0, 1, "S-0-0386") == 'ns=7;s=Sercos,0,1.ParameterSet."S-0-0386"'


def test_plugin(server):
    with MonitoringSource.from_name("OPCUAMonitoringSource", subarray=SUBARRAY, **pointing_config(server)) as source:
        assert isinstance(source, OPCUAMonitoringSource)
        assert source.url == f"opc.tcp://127.0.0.1:{server.port}"
        assert source.monitoring_types == (MonitoringType.TELESCOPE_POINTINGS,)


def test_fill_monitoring_container(server):
    server.set(1, 70.0)
    server.set(2, 180.0)
    with OPCUAMonitoringSource(subarray=SUBARRAY, **pointing_config(server)) as source:
        assert source.is_connected
        # event time in TAI, as the events of the camera
        event = event_at(Time.now().tai)
        source.fill_monitoring_container(event)
        pointing = event.monitoring.tel[TEL_ID].pointing
        assert u.isclose(pointing.altitude, 70 * u.deg)
        assert u.isclose(pointing.azimuth, 180 * u.deg)

        # the values follow the server
        server.set(1, 45.5)
        assert wait_for(lambda: u.isclose(source.values_at(Time.now()).get("pointing.altitude", np.nan * u.deg), 45.5 * u.deg))
        values = source.get_values(MonitoringType.TELESCOPE_POINTINGS, Time.now(), tel_id=TEL_ID)
        assert u.isclose(values["altitude"], 45.5 * u.deg)
        # the past events get the values read at their time
        event = event_at(Time(event.trigger.time))
        source.fill_monitoring_container(event)
        assert u.isclose(event.monitoring.tel[TEL_ID].pointing.altitude, 70 * u.deg)

        # values of several times, NaN out of the time range of the buffer (timestamp_tolerance)
        times = Time.now() + [-1, 0, 1] * u.hour
        values = source.get_values(MonitoringType.TELESCOPE_POINTINGS, times)
        assert np.isnan(values["altitude"][[0, 2]]).all()
        assert u.isclose(values["altitude"][1], 45.5 * u.deg)
        event = event_at(times[0])
        source.fill_monitoring_container(event)
        assert np.isnan(event.monitoring.tel[TEL_ID].pointing.altitude)

        table = source.get_table(MonitoringType.TELESCOPE_POINTINGS, tel_id=TEL_ID)
        assert len(table) > 1
        assert table["pointing.altitude"].unit == u.deg
        assert np.all(np.diff(table["time"].unix) > 0)
        with pytest.raises(KeyError):
            source.get_values(MonitoringType.CAMERA_COEFFICIENTS, Time.now())
        with pytest.raises(KeyError):
            source.get_table(MonitoringType.TELESCOPE_POINTINGS, tel_id=22)


def test_scale_offset_and_unreadable_node(server):
    server.set(1, 1000.0)
    config = pointing_config(server)
    config["nodes"]["pointing.azimuth"] = server.node_id(2, "S-0-9999")  # unknown node
    with OPCUAMonitoringSource(
        subarray=SUBARRAY, scales={"pointing.altitude": 1e-2}, offsets={"pointing.altitude": 5.0}, **config,
    ) as source:
        event = event_at(Time.now())
        source.fill_monitoring_container(event)
        pointing = event.monitoring.tel[TEL_ID].pointing
        assert u.isclose(pointing.altitude, 15 * u.deg)
        # not filled
        assert np.isnan(pointing.azimuth)


def test_unusable_value_and_clock_offset(server):
    server.set(2, 180.0)
    config = pointing_config(server)
    # text value of a node with a unit: only this field is missing
    config["nodes"]["pointing.altitude"] = f"ns={server.namespace};s=status"
    with OPCUAMonitoringSource(subarray=SUBARRAY, **config) as source:
        assert source.is_connected
        event = event_at(Time.now())
        source.fill_monitoring_container(event)
        pointing = event.monitoring.tel[TEL_ID].pointing
        assert np.isnan(pointing.altitude)
        assert u.isclose(pointing.azimuth, 180 * u.deg)
        assert source.is_connected
        assert source.n_events_without_values == 0
        # events whose clock is 1 min apart from the computer: counted
        source.fill_monitoring_container(event_at(Time.now() - 1 * u.min))
        assert source.n_events_without_values == 1


def test_buffer_order_and_arrays(server):
    with OPCUAMonitoringSource(
        subarray=SUBARRAY, host="127.0.0.1", port=server.port, polling_period=60, buffer_duration=10,
        nodes={"pointing.altitude": server.node_id(1)},
    ) as source:
        source._times.clear()
        source._samples.clear()
        now = Time.now().unix
        # clock stepped back: kept in time order; older than buffer_duration: removed
        for offset, value in ((0, 1.0), (-5, 2.0), (2, 3.0), (-20, 4.0)):
            source._add_sample(now + offset, [value])
        assert source._times == [now - 5, now, now + 2]
        times = Time([now - 5, now, now + 2, now + 100], format="unix")
        values = source.get_values(MonitoringType.TELESCOPE_POINTINGS, times)
        # without unit: numpy array
        assert isinstance(values["altitude"], np.ndarray)
        np.testing.assert_array_equal(values["altitude"], [2.0, 1.0, 3.0, np.nan])


def test_server_not_reachable():
    with OPCUAMonitoringSource(
        subarray=SUBARRAY, host="127.0.0.1", port=free_port(), polling_period=POLLING_PERIOD, timeout=0.5,
        nodes={"pointing.altitude": "ns=2;s=altitude"}, units={"pointing.altitude": "deg"},
    ) as source:
        assert not source.is_connected
        event = event_at(Time.now())
        source.fill_monitoring_container(event)
        assert np.isnan(event.monitoring.tel[TEL_ID].pointing.altitude)


def test_reconnection():
    port = free_port()
    drive = DriveServer(port)
    config = pointing_config(drive, reconnect_period=0.2, timeout=0.5)
    try:
        with OPCUAMonitoringSource(subarray=SUBARRAY, **config) as source:
            assert source.is_connected
            drive.stop()
            assert wait_for(lambda: not source.is_connected)
            drive = DriveServer(port)
            drive.set(1, 12.0)
            assert wait_for(lambda: source.values_at(Time.now()).get("pointing.altitude") == 12 * u.deg, timeout=10)
    finally:
        drive.stop()


def test_configuration_errors(server):
    with pytest.raises(ValueError, match="not a field"):
        OPCUAMonitoringSource(subarray=SUBARRAY, nodes={"pointing.elevation": "ns=2;s=x"})
    with pytest.raises(ValueError, match="without node"):
        OPCUAMonitoringSource(subarray=SUBARRAY, nodes={}, units={"pointing.altitude": "deg"})
    with pytest.raises(ValueError, match="tel_id"):
        OPCUAMonitoringSource(subarray=SUBARRAY_DESCRIPTION, nodes={})
    with pytest.raises(ValueError, match="not in the subarray"):
        OPCUAMonitoringSource(subarray=SUBARRAY, tel_id=22, nodes={})
    with pytest.raises(ValueError, match="is a container"):
        OPCUAMonitoringSource(subarray=SUBARRAY, nodes={"pointing": "ns=2;s=x"})
