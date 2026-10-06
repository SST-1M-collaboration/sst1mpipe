"""
Files distributed with sst1mpipe:

- ``config``: configuration files of the analysis (data, MC at low and high NSB, RTA)
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
    "MC_CONFIG_FILES",
    "RTA_CONFIG_FILE",
    "PDE_CORRECTION_FACTORS_FILE",
    "SUBARRAY_FILE",
    "PIXEL_MAPPING_FILE",
    "CAMERA_CONFIG_FILE",
]

RESOURCES_DIR = files(__name__)

CONFIG_DIR = RESOURCES_DIR / "config"
CALIBRATION_DIR = RESOURCES_DIR / "calibration"
WINDOW_DIR = CALIBRATION_DIR / "window"
DARK_RUNS_DIR = CALIBRATION_DIR / "dark_runs"
INSTRUMENT_DIR = RESOURCES_DIR / "instrument"
TEST_DATA_DIR = RESOURCES_DIR / "test_data"

DATA_CONFIG_FILE = CONFIG_DIR / "sst1mpipe_data_config.json"
MC_CONFIG_FILES = {nsb: CONFIG_DIR / f"sst1mpipe_mc_config_{nsb}_nsb.json" for nsb in ("low", "high")}
RTA_CONFIG_FILE = CONFIG_DIR / "sst1mpipe_rta_config.json"

PDE_CORRECTION_FACTORS_FILE = CALIBRATION_DIR / "mc_pde_correction_factors.json"

SUBARRAY_FILE = INSTRUMENT_DIR / "sst1m_array.h5"
PIXEL_MAPPING_FILE = INSTRUMENT_DIR / "digicam_pixels_mapping_V5T.txt"
CAMERA_CONFIG_FILE = INSTRUMENT_DIR / "camera_config.cfg"
