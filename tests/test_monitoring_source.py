"""
SST1MMonitoringSource (monitoring files of sst1mpipe-monitoring) and sst1mpipe-process with
a monitoring source
"""
import astropy.units as u
import numpy as np
import pytest
import tables
from ctapipe.core import run_tool
from ctapipe.io import MonitoringSource, read_table
from ctapipe.io.hdf5dataformat import DL1_PIXEL_STATISTICS_GROUP
from ctapipe.io.monitoringtypes import MonitoringType

from sst1mpipe.io import R0_PIXEL_STATISTICS_GROUP, get_dl1_info, load_dl1_sst1m
from sst1mpipe.io.containers import PEDESTAL_EVENT_TYPES
from sst1mpipe.io.monitoring_source import SST1MMonitoringSource, statistics_container
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.resources import MONITORING_CONFIG_FILE, RTA_CONFIG_FILE, TEST_DATA_DIR
from sst1mpipe.scripts.sst1mpipe_monitoring_tool import MonitoringTool
from sst1mpipe.scripts.sst1mpipe_process_tool import ProcessorTool

# observation of Mrk 421 by tel 22: 35 events, 13 pedestal events
INPUT_FILE = TEST_DATA_DIR / "zfits" / "SST1M2_20260120_1102.fits.fz"
TEL_ID = 22
N_EVENTS, N_PEDESTALS = 35, 13


def monitoring_file(path, *args):
    run_tool(MonitoringTool(), argv=[
        f"--input={INPUT_FILE}",
        f"--output={path}",
        f"--config={MONITORING_CONFIG_FILE}",
        "--no-progress",
        *args,
    ], raises=True)
    return path


def process(input_file, output, *args):
    tool = ProcessorTool()
    run_tool(tool, argv=[
        f"--input={input_file}",
        f"--output={output}",
        f"--config={RTA_CONFIG_FILE}",
        "--no-progress",
        "--ProcessorTool.wobble_in_output_name=False",
        *args,
    ], raises=True)
    return tool


@pytest.fixture(scope="module")
def monitoring(tmp_path_factory):
    path = tmp_path_factory.mktemp("monitoring") / "run.monitoring.h5"
    return monitoring_file(path, "--WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_size=5")


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
            source_pointing = event.monitoring.tel[TEL_ID].pointing
            source_altitude, source_azimuth = source_pointing.altitude, source_pointing.azimuth
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
            assert u.isclose(monitoring.pointing.altitude, source_altitude, atol=0.01 * u.deg)
            assert u.isclose(monitoring.pointing.azimuth, source_azimuth, atol=0.01 * u.deg)


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
    """the same events and images as with the pedestal monitors of the processing"""
    output = tmp_path / "events.dl1.h5"
    tool = process(
        INPUT_FILE, output,
        "--monitoring-source=SST1MMonitoringSource",
        f"--SST1MMonitoringSource.input_files={monitoring}",
    )
    assert tool.r0_pedestal_monitor is None and tool.dl1_pedestal_monitor is None
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
    data = load_dl1_sst1m(str(output), tel=f"tel_{TEL_ID:03d}")
    assert np.all((70 < data["true_alt_tel"]) & (data["true_alt_tel"] < 80))


def test_process_without_attached_monitoring(monitoring, tmp_path):
    output = tmp_path / "events.dl1.h5"
    process(
        INPUT_FILE, output,
        "--monitoring-source=SST1MMonitoringSource",
        f"--SST1MMonitoringSource.input_files={monitoring}",
        "--ProcessorTool.attach_monitoring=False",
    )
    with tables.open_file(output) as dl1:
        assert "/dl0/monitoring" not in dl1 and "/r0/monitoring" not in dl1 and "/dl1/monitoring" not in dl1


def test_events_without_monitoring(tmp_path):
    """
    monitoring of the 20 first events, without timestamp tolerance: the events out of the chunks
    of pedestal events or of the pointing table are skipped
    """
    monitoring = monitoring_file(tmp_path / "run.monitoring.h5", "--max-events=20")
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
    tool = process(
        INPUT_FILE, output,
        "--monitoring-source=SST1MMonitoringSource",
        f"--SST1MMonitoringSource.input_files={monitoring}",
        "--SST1MMonitoringSource.timestamp_tolerance=0 s",
    )
    assert tool.n_events_without_monitoring == n_without
    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{TEL_ID:03d}")
    assert len(parameters) == N_EVENTS - N_PEDESTALS - n_without
    assert get_dl1_info(output)[f"n_triggered_tel{TEL_ID - 20}"][0] == N_EVENTS
