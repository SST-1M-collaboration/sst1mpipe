from importlib.resources import files

import numpy as np
import pytest
from ctapipe.core import run_tool
from ctapipe.image.cleaning import NSBImageCleaner
from ctapipe.io import read_table

from sst1mpipe.calib.calib import DEFAULT_CALIBRATION_FILES
from sst1mpipe.io import get_dl1_info
from sst1mpipe.scripts.sst1mpipe_process_tool import ProcessorTool

FILE_TEL_1 = files('sst1mpipe.resources.zfits').joinpath('SST1M1_20260121_0001.fits.fz')
RTA_CONFIG = files('sst1mpipe.data').joinpath('sst1mpipe_rta_config.json')
N_EVENTS = 150


@pytest.mark.parametrize("voltage_drop_correction", ["global", "none"])
def test_process_r0_file(tmp_path, voltage_drop_correction):

    output = tmp_path / "events.dl1.h5"
    tool = ProcessorTool()
    ret = run_tool(tool, argv=[
        f"--input={FILE_TEL_1}",
        f"--output={output}",
        f"--config={RTA_CONFIG}",
        f"--max-events={N_EVENTS}",
        "--ProcessorTool.progress_bar=False",
        f"--R0R1Calibrator.voltage_drop_correction={voltage_drop_correction}",
        "--R0PedestalMonitor.n_events=10",
    ], raises=True)

    assert ret == 0
    # the R0 -> R1 calibration is configured from the config file and the command line
    assert tool.r0_r1_calibrator.voltage_drop_correction.tel[21] == voltage_drop_correction
    assert tool.r0_pedestal_monitor.n_events.tel[21] == 10
    assert tool.r0_pedestal_monitor.n_buffered(21) == 10

    parameters = read_table(output, "/dl1/event/telescope/parameters/tel_021")
    assert len(parameters) == N_EVENTS
    assert np.isfinite(parameters["camera_frame_hillas_intensity"]).sum() == 0  # dark run

    info = get_dl1_info(output)
    assert info["calib_file"][0].endswith(DEFAULT_CALIBRATION_FILES[21])
    assert info["n_pedestal"][0] == N_EVENTS
    assert info["n_triggered_tel1"][0] == N_EVENTS


def test_dl1_pedestal_monitor_used_by_nsb_image_cleaner(tmp_path, monkeypatch):

    # record the pedestal std given to the NSBImageCleaner for each event
    pedestal_stds = []
    nsb_image_cleaner_call = NSBImageCleaner.__call__

    def record(self, tel_id, image, arrival_times=None, *, monitoring=None):
        std = monitoring.pedestal.charge_std
        pedestal_stds.append(None if std is None else np.array(std))
        return nsb_image_cleaner_call(self, tel_id, image, arrival_times, monitoring=monitoring)

    monkeypatch.setattr(NSBImageCleaner, "__call__", record)

    tool = ProcessorTool()
    run_tool(tool, argv=[
        f"--input={FILE_TEL_1}",
        f"--output={tmp_path / 'events.dl1.h5'}",
        f"--config={RTA_CONFIG}",
        f"--max-events={N_EVENTS}",
        "--ProcessorTool.progress_bar=False",
        "--ImageProcessor.image_cleaner_type=NSBImageCleaner",
        "--DL1PedestalMonitor.n_events=20",
        # in the dark run of the test file all the pixels are dead (pedestal std < 2.5 ADC):
        # their image is 0, so they are kept to have a pedestal std
        "--R0R1Calibrator.flag_dead_pixels=False",
    ], raises=True)

    # all the events of the test file are pedestal events
    monitor = tool.dl1_pedestal_monitor
    assert monitor.processed_events[21] == N_EVENTS
    assert monitor.n_buffered(21) == 20

    # no statistics before the first pedestal event, then the std of the images of the sliding window
    assert len(pedestal_stds) == N_EVENTS
    assert pedestal_stds[0] is None
    assert all(std is not None and std.shape == (1296, ) for std in pedestal_stds[1:])
    assert np.all(np.isfinite(pedestal_stds[-1]))
    # the std of the images (p.e.) grows from 0 (a single image in the window)
    assert np.all(pedestal_stds[1] == 0)
    assert 0.1 < np.median(pedestal_stds[-1]) < 2
