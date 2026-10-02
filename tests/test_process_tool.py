from importlib.resources import files

import numpy as np
import pytest
from ctapipe.core import run_tool
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
