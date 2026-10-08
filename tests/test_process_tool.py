import numpy as np
import pytest
from ctapipe.core import run_tool
from ctapipe.image.cleaning import NSBImageCleaner
from ctapipe.io import read_table

from sst1mpipe.calib.calib import DEFAULT_CALIBRATION_FILES
from sst1mpipe.io import get_dl1_info
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

    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
    assert len(parameters) == run["n_events"]
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
        std = monitoring.pedestal.charge_std
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
