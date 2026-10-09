import astropy.units as u
import numpy as np
import pytest
from ctapipe.core import run_tool
from ctapipe.image.cleaning import NSBImageCleaner
from ctapipe.io import HDF5EventSource, read_table
from ctapipe.monitoring import PedestalImageInterpolator

from sst1mpipe.calib.calib import DEFAULT_CALIBRATION_FILES
from sst1mpipe.io import (
    DL1_PEDESTAL_GROUP,
    R0_PEDESTAL_GROUP,
    get_dl1_info,
    load_dl1_pedestals,
    load_dl1_sst1m,
    load_r0_pedestals,
)
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.scripts.sst1mpipe_process_tool import ProcessorTool
from sst1mpipe.resources import RTA_CONFIG_FILE, TEST_DATA_DIR

# observation of Mrk 421 (wobble W1) by both telescopes, with Cherenkov and pedestal events
OBSERVATION_FILES = {
    21: dict(path=TEST_DATA_DIR / "zfits" / "SST1M1_20260120_1179.fits.fz", n_events=150, n_pedestals=40),
    22: dict(path=TEST_DATA_DIR / "zfits" / "SST1M2_20260120_1102.fits.fz", n_events=35, n_pedestals=13),
}
TARGET, WOBBLE, RA, DEC = "MRK421", "W1", 166.9948, 38.1053


@pytest.mark.parametrize("tel_id", sorted(OBSERVATION_FILES))
@pytest.mark.parametrize("voltage_drop_correction", ["global", "none"])
def test_process_r0_file(tmp_path, tel_id, voltage_drop_correction):
    """regular processing of an observation run, from R0 up to the DL1 parameters"""
    run = OBSERVATION_FILES[tel_id]
    tool = ProcessorTool()
    ret = run_tool(tool, argv=[
        f"--input={run['path']}",
        f"--output={tmp_path / 'events.dl1.h5'}",
        f"--config={RTA_CONFIG_FILE}",
        f"--max-events={run['n_events']}",
        "--ProcessorTool.progress_bar=False",
        f"--R0R1Calibrator.voltage_drop_correction={voltage_drop_correction}",
        "--R0PedestalMonitor.n_events=10",
    ], raises=True)

    assert ret == 0
    # the R0 -> R1 calibration is configured from the config file and the command line
    assert tool.r0_r1_calibrator.voltage_drop_correction.tel[tel_id] == voltage_drop_correction
    assert tool.r0_pedestal_monitor.n_buffered(tel_id) == 10
    assert tool.n_skipped_events == 0

    # the wobble of the run is added to the name of the output file
    output = tmp_path / f"events_{WOBBLE}.dl1.h5"
    assert output.exists()

    # the pedestal events are not written as events
    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
    assert len(parameters) == run["n_events"] - run["n_pedestals"]
    intensity = parameters["camera_frame_hillas_intensity"]
    assert np.isfinite(intensity).sum() > 0
    assert np.all(intensity[np.isfinite(intensity)] > 0)

    info = get_dl1_info(output)
    # calibration file of the telescope
    assert info["calib_file"][0].endswith(DEFAULT_CALIBRATION_FILES[tel_id])
    assert info[f"n_triggered_tel{tel_id - 20}"][0] == run["n_events"]
    assert info["n_pedestal"][0] == run["n_pedestals"]
    # the pedestal events do not survive the cleaning
    assert info["n_survived_pedestals"][0] == 0
    assert (info["target"][0], info["wobble"][0]) == (TARGET, WOBBLE)
    assert (info["ra"][0], info["dec"][0]) == pytest.approx((RA, DEC))
    assert not info["manual_coords"][0]


def test_dl1_pedestal_monitor_used_by_nsb_image_cleaner(tmp_path, monkeypatch):

    # record the pedestal std given to the NSBImageCleaner for each event
    pedestal_stds = []
    nsb_image_cleaner_call = NSBImageCleaner.__call__

    def record(self, tel_id, image, arrival_times=None, *, monitoring=None):
        std = monitoring.pixel_statistics.pedestal_image.std
        pedestal_stds.append(None if std is None else np.array(std))
        return nsb_image_cleaner_call(self, tel_id, image, arrival_times, monitoring=monitoring)

    monkeypatch.setattr(NSBImageCleaner, "__call__", record)

    run = OBSERVATION_FILES[21]
    tool = ProcessorTool()
    run_tool(tool, argv=[
        f"--input={run['path']}",
        f"--output={tmp_path / 'events.dl1.h5'}",
        f"--config={RTA_CONFIG_FILE}",
        f"--max-events={run['n_events']}",
        "--ProcessorTool.progress_bar=False",
        "--DL1PedestalMonitor.n_events=20",
    ], raises=True)

    monitor = tool.dl1_pedestal_monitor
    assert monitor.processed_events[21] == run["n_pedestals"]
    assert monitor.n_buffered(21) == 20

    # no statistics before the first pedestal event (the second event of the run),
    # then the std of the images of the sliding window
    assert len(pedestal_stds) == run["n_events"]
    assert pedestal_stds[0] is None and pedestal_stds[1] is None
    assert all(std is not None and std.shape == (1296, ) for std in pedestal_stds[2:])
    assert np.all(np.isfinite(pedestal_stds[-1]))
    # the std of the images (p.e.) grows from 0 (a single image in the window)
    assert np.all(pedestal_stds[2] == 0)
    assert 0.1 < np.median(pedestal_stds[-1]) < 2


def test_monitoring_tables(tmp_path):
    """pointing and pedestal monitoring tables written by sst1mpipe-process"""
    tel_id, run = 22, OBSERVATION_FILES[22]
    tel = f"tel_{tel_id:03d}"
    output = tmp_path / "events.dl1.h5"
    run_tool(ProcessorTool(), argv=[
        f"--input={run['path']}",
        f"--output={output}",
        f"--config={RTA_CONFIG_FILE}",
        "--ProcessorTool.wobble_in_output_name=False",
        "--ProcessorTool.pedestal_monitoring_interval=5",
        "--R0PedestalMonitor.n_events=4",
    ], raises=True)

    # pointing of the events, from the event source
    with SST1MEventSource(run["path"]) as source:
        pointing = {
            event.index.event_id: (event.trigger.time, event.monitoring.tel[tel_id].pointing.azimuth,
                                   event.monitoring.tel[tel_id].pointing.altitude)
            for event in source
        }
    times = [time for time, _, _ in pointing.values()]

    # a row each time the alt/az is computed (every second) and at the time of the last event
    table = read_table(output, f"/dl0/monitoring/telescope/pointing/{tel}")
    assert table.colnames == ["time", "azimuth", "altitude"]
    assert np.all(np.diff(table["time"].mjd) > 0)
    assert table["time"][0] == times[0]
    assert table["time"][-1] == times[-1]
    assert 1 < len(table) < run["n_events"]

    # the pointing of the events is interpolated by ctapipe (HDF5EventSource, load_dl1_sst1m)
    with HDF5EventSource(output) as source:
        interpolated = {e.index.event_id: e.monitoring.tel[tel_id].pointing for e in source}
    assert len(interpolated) == run["n_events"] - run["n_pedestals"]
    for event_id, p in interpolated.items():
        _, azimuth, altitude = pointing[event_id]
        assert u.isclose(p.azimuth, azimuth, atol=0.01 * u.deg)
        assert u.isclose(p.altitude, altitude, atol=0.01 * u.deg)
    data = load_dl1_sst1m(str(output), tel=tel)
    np.testing.assert_allclose(data["true_az_tel"], [interpolated[e].azimuth.to_value(u.deg) for e in data["event_id"]])
    np.testing.assert_allclose(data["true_alt_tel"], [interpolated[e].altitude.to_value(u.deg) for e in data["event_id"]])

    # statistics of the pedestal events, every 5 pedestal events, in chunks of pedestal events
    n_rows = run["n_pedestals"] // 5
    r0 = read_table(output, f"{R0_PEDESTAL_GROUP}/{tel}")
    dl1 = read_table(output, f"{DL1_PEDESTAL_GROUP}/{tel}")
    assert DL1_PEDESTAL_GROUP == "/dl1/monitoring/telescope/calibration/camera/pixel_statistics/sky_pedestal_image"
    assert len(r0) == len(dl1) == n_rows
    # last pedestal event of the chunks: the same for both
    np.testing.assert_array_equal(r0["event_id_end"], dl1["event_id_end"])
    assert np.all(r0["time_start"] <= r0["time_end"])
    # sliding windows of 4 (R0PedestalMonitor) and 1000 (DL1PedestalMonitor) pedestal events
    assert list(r0["n_events"]) == [4] * n_rows
    assert list(dl1["n_events"]) == [5 * (i + 1) for i in range(n_rows)]
    assert np.all(dl1["time_start"] == dl1["time_start"][0])
    assert r0["std"].shape == (n_rows, 1, 1296)
    assert dl1["std"].shape == (n_rows, 1296)
    assert np.all(dl1["is_valid"]) and not np.any(dl1["outlier_mask"])
    # ADC baseline (~ 250-350 ADC) and std of the images (~ 1 p.e.)
    assert 100 < np.nanmedian(r0["mean"]) < 1000
    assert 0.1 < np.nanmedian(dl1["std"]) < 3
    # same tables with the sst1mpipe readers
    assert len(load_r0_pedestals(output)) == len(load_dl1_pedestals(output, tel=tel)) == n_rows

    # the statistics of the images are interpolated by ctapipe, as in the HDF5MonitoringSource
    interpolator = PedestalImageInterpolator()
    interpolator.add_table(tel_id, dl1)
    np.testing.assert_allclose(interpolator(tel_id, dl1["time_end"][-1])["std"], dl1["std"][-1])
