"""
Files distributed with sst1mpipe:

- ``config``: configuration files of the analysis (data, MC, RTA, trigger emulator)
- ``calibration``: camera calibration (dc_to_pe from dark runs, window transmittance,
  PDE drop correction of the simulations, fits of dark runs)
- ``instrument``: description of the telescopes and of the camera (subarray, pixel mapping)
- ``test_data``: raw data files used by the tests (git LFS)
"""
from importlib.resources import files

__all__ = [
    "RESOURCES_DIR",
    "CONFIG_DIR",
    "CALIBRATION_DIR",
    "WINDOW_DIR",
    "DARK_RUNS_DIR",
    "INSTRUMENT_DIR",
    "TEST_DATA_DIR",
    "DATA_CONFIG_FILE",
    "MC_CONFIG_FILE",
    "RTA_CONFIG_FILE",
    "PDE_CORRECTION_FACTORS_FILE",
    "SUBARRAY_FILE",
    "PIXEL_MAPPING_FILE",
    "CAMERA_CONFIG_FILE",
    "SST1M_POINTING_MODEL_CONFIG_FILE",
]

RESOURCES_DIR = files(__name__)

CONFIG_DIR = RESOURCES_DIR / "config"
CALIBRATION_DIR = RESOURCES_DIR / "calibration"
WINDOW_DIR = CALIBRATION_DIR / "window"
DARK_RUNS_DIR = CALIBRATION_DIR / "dark_runs"
INSTRUMENT_DIR = RESOURCES_DIR / "instrument"
TEST_DATA_DIR = RESOURCES_DIR / "test_data"

DATA_CONFIG_FILE = CONFIG_DIR / "sst1mpipe_data_config.json"
MC_CONFIG_FILE = CONFIG_DIR / "sst1mpipe_mc_config.json"
RTA_CONFIG_FILE = CONFIG_DIR / "sst1mpipe_rta_config.json"
# TPointPointingModel of the SST-1M drives (AzElCorrCoefTel<n>.conf), periods from noon to noon UTC: tel 21 from
# the night 20230420 (same coefficients in its three files, PMAz2 and PMEl2 from 0 to -1 on 20230612); tel 22 for
# the nights 20240430 to 20260313 and from 20260317 (new mirrors). The other times use the latest coefficients of
# the telescope, with a warning.
SST1M_POINTING_MODEL_CONFIG_FILE = CONFIG_DIR / "sst1m_tpoint_pointing_model.json"

PDE_CORRECTION_FACTORS_FILE = CALIBRATION_DIR / "mc_pde_correction_factors.json"

SUBARRAY_FILE = INSTRUMENT_DIR / "sst1m_array.h5"
PIXEL_MAPPING_FILE = INSTRUMENT_DIR / "digicam_pixels_mapping_V5T.txt"
CAMERA_CONFIG_FILE = INSTRUMENT_DIR / "camera_config.cfg"
