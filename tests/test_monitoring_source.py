"""
SST1MMonitoringSource (monitoring files of sst1mpipe-monitoring) and sst1mpipe-process with
a monitoring source
"""
import astropy.units as u
import numpy as np
import pytest
import tables
from astropy.coordinates import AltAz
from ctapipe.io import HDF5EventSource, MonitoringSource, read_table
from ctapipe.io.hdf5dataformat import DL1_PIXEL_STATISTICS_GROUP
from ctapipe.io.monitoringtypes import MonitoringType

from conftest import OBSERVATION_RUNS, monitoring_source_args, run_monitoring, run_process
from sst1mpipe.io import R0_PIXEL_STATISTICS_GROUP, get_dl1_info, load_dl1_sst1m
from sst1mpipe.io.containers import PEDESTAL_EVENT_TYPES
from sst1mpipe.io.monitoring_source import SST1MMonitoringSource, statistics_container
from sst1mpipe.io.sst1m_event_source import SST1MEventSource

# observation of Mrk 421 by tel 22
TEL_ID = 22
INPUT_FILE = OBSERVATION_RUNS[TEL_ID]["path"]
N_EVENTS, N_PEDESTALS = OBSERVATION_RUNS[TEL_ID]["n_events"], OBSERVATION_RUNS[TEL_ID]["n_pedestals"]


@pytest.fixture(scope="module")
def monitoring(observation_monitoring):
    """monitoring file of the observation of tel 22, chunks of 5 pedestal events"""
    return observation_monitoring[TEL_ID]


def pointing_altaz(source, time):
    """alt/az of the pointing of the run (RA/Dec of the event source) at the site of the array"""
    return source.pointing.transform_to(AltAz(obstime=time, location=source.subarray.reference_location))


def test_monitoring_source_plugin(monitoring):
    with SST1MEventSource(INPUT_FILE) as source:
        monitoring_source = MonitoringSource.from_name(
            "SST1MMonitoringSource", subarray=source.subarray, input_files=[monitoring],
        )
        # the telescope of the run
        assert source.subarray.tel_ids.tolist() == monitoring_source.subarray.tel_ids.tolist() == [TEL_ID]
        assert isinstance(monitoring_source, SST1MMonitoringSource)


def test_fill_monitoring_container(monitoring):
    waveform = read_table(monitoring, f"{R0_PIXEL_STATISTICS_GROUP}/sky_pedestal_waveform/tel_{TEL_ID:03d}")
    dead_pixels = np.all(waveform["outlier_mask"][:, 0], axis=0)
    assert dead_pixels.sum() > 0
    with SST1MEventSource(INPUT_FILE) as source:
        monitoring_source = SST1MMonitoringSource(subarray=source.subarray, input_files=[monitoring])
        for event in source:
            # the event source gives the RA/Dec of the pointing only
            assert np.isnan(event.monitoring.tel[TEL_ID].pointing.altitude)
            altaz = pointing_altaz(source, event.trigger.time)
            monitoring_source.fill_monitoring_container(event)
            monitoring = event.monitoring.tel[TEL_ID]
            statistics = monitoring.camera.pixel_statistics

            # statistics of the ADC samples, (n_channels, n_pixels): the outliers of all the chunks
            # (dead pixels) are NaN
            std = statistics.pedestal_waveform.std
            assert std.shape == (1, 1296)
            assert np.all(np.isnan(std[0, dead_pixels]))
            np.testing.assert_allclose(np.nanmedian(std), np.nanmedian(waveform["std"]), rtol=0.05)
            # statistics of the images, channel 0
            assert statistics.pedestal_image.std.shape == (1296, )
            # pointing interpolated from the pointing table
            assert u.isclose(monitoring.pointing.altitude, altaz.alt, atol=0.01 * u.deg)
            assert u.isclose(monitoring.pointing.azimuth, altaz.az, atol=0.01 * u.deg)


def test_no_statistics_out_of_the_chunks(monitoring):
    """no statistics without chunk of pedestal events at the time of the event (timestamp_tolerance)"""
    monitoring_source = SST1MMonitoringSource(input_files=[monitoring], timestamp_tolerance=0 * u.s)
    image = read_table(monitoring, f"{DL1_PIXEL_STATISTICS_GROUP}/sky_pedestal_image/tel_{TEL_ID:03d}")

    def statistics(time):
        values = monitoring_source.get_values(
            MonitoringType.PIXEL_STATISTICS, time=time, tel_id=TEL_ID, subtype="sky_pedestal_image",
        )
        return statistics_container(values, channel=0)

    assert statistics(image["time_start"][0]).std.shape == (1296, )
    assert statistics(image["time_start"][0] - 1 * u.s).std is None
    assert statistics(image["time_end"][-1] + 1 * u.s).std is None


def test_process_with_monitoring_source(monitoring, tmp_path):
    """the monitoring data are used for the calibration and attached to the DL1 file"""
    output = tmp_path / "events.dl1.h5"
    tool = run_process(INPUT_FILE, output, *monitoring_source_args(monitoring))
    assert tool.n_events_without_monitoring == 0

    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{TEL_ID:03d}")
    assert len(parameters) == N_EVENTS - N_PEDESTALS
    intensity = parameters["camera_frame_hillas_intensity"]
    assert np.isfinite(intensity).sum() > 0

    info = get_dl1_info(output)
    assert info["n_pedestal"][0] == N_PEDESTALS
    assert info[f"n_triggered_tel{TEL_ID - 20}"][0] == N_EVENTS

    # the monitoring data are copied into the DL1 file, the pointing is read by load_dl1_sst1m
    with tables.open_file(monitoring) as source, tables.open_file(output) as dl1:
        monitoring_tables = [t._v_pathname for t in source.walk_nodes("/", "Table") if "monitoring" in t._v_pathname]
        assert len(monitoring_tables) == 3
        for path in monitoring_tables:
            np.testing.assert_array_equal(dl1.get_node(path).read(), source.get_node(path).read())
    # the pointing of the events is interpolated by ctapipe (HDF5EventSource) and by load_dl1_sst1m
    with SST1MEventSource(INPUT_FILE) as source:
        times = {e.index.event_id: e.trigger.time for e in source}
        altaz = {event_id: pointing_altaz(source, time) for event_id, time in times.items()}
    with HDF5EventSource(output) as source:
        interpolated = {e.index.event_id: e.monitoring.tel[TEL_ID].pointing for e in source}
    assert len(interpolated) == N_EVENTS - N_PEDESTALS
    for event_id, pointing in interpolated.items():
        assert u.isclose(pointing.altitude, altaz[event_id].alt, atol=0.01 * u.deg)
        assert u.isclose(pointing.azimuth, altaz[event_id].az, atol=0.01 * u.deg)
    data = load_dl1_sst1m(str(output), tel=f"tel_{TEL_ID:03d}")
    np.testing.assert_allclose(data["true_alt_tel"], [interpolated[e].altitude.to_value(u.deg) for e in data["event_id"]])
    np.testing.assert_allclose(data["true_az_tel"], [interpolated[e].azimuth.to_value(u.deg) for e in data["event_id"]])


def test_process_without_attached_monitoring(monitoring, tmp_path):
    output = tmp_path / "events.dl1.h5"
    run_process(INPUT_FILE, output, *monitoring_source_args(monitoring), "--ProcessorTool.attach_monitoring=False")
    with tables.open_file(output) as dl1:
        assert "/dl0/monitoring" not in dl1 and "/r0/monitoring" not in dl1 and "/dl1/monitoring" not in dl1


def test_events_without_monitoring(tmp_path):
    """
    monitoring of the 20 first events, without timestamp tolerance: the events out of the chunks
    of pedestal events or of the pointing table are skipped
    """
    monitoring = tmp_path / "run.monitoring.h5"
    run_monitoring(INPUT_FILE, monitoring, "--max-events=20")
    chunks = [
        read_table(monitoring, f"{group}/tel_{TEL_ID:03d}")
        for group in (f"{R0_PIXEL_STATISTICS_GROUP}/sky_pedestal_waveform", f"{DL1_PIXEL_STATISTICS_GROUP}/sky_pedestal_image")
    ]
    pointing = read_table(monitoring, f"/dl0/monitoring/telescope/pointing/tel_{TEL_ID:03d}")

    def has_monitoring(time):
        in_chunks = all(np.any((table["time_start"] <= time) & (time <= table["time_end"])) for table in chunks)
        return in_chunks and pointing["time"][0] <= time <= pointing["time"][-1]

    with SST1MEventSource(INPUT_FILE) as source:
        events = [(e.trigger.time, e.trigger.event_type) for e in source]
    n_without = sum(not has_monitoring(time) for time, event_type in events if event_type not in PEDESTAL_EVENT_TYPES)
    assert n_without > 0

    output = tmp_path / "events.dl1.h5"
    tool = run_process(
        INPUT_FILE, output, *monitoring_source_args(monitoring), "--SST1MMonitoringSource.timestamp_tolerance=0 s",
    )
    assert tool.n_events_without_monitoring == n_without
    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{TEL_ID:03d}")
    assert len(parameters) == N_EVENTS - N_PEDESTALS - n_without
    assert get_dl1_info(output)[f"n_triggered_tel{TEL_ID - 20}"][0] == N_EVENTS
