"""
Parameters of the analysis computed when the DL1 files of sst1mpipe-process are read
(load_dl1_sst1m, add_dl1_extra_parameters)
"""
from copy import deepcopy

import astropy.units as u
import numpy as np
import pytest
from ctapipe.instrument import SubarrayDescription
from ctapipe.io import read_table
from ctapipe.io.hdf5dataformat import DL1_PIXEL_STATISTICS_GROUP

from conftest import OBSERVATION_RUNS, monitoring_source_args, run_process
from sst1mpipe.io import load_config, load_dl1_sst1m, read_trigger_time
from sst1mpipe.resources import DATA_CONFIG_FILE, MC_CONFIG_FILES, TEST_DATA_DIR

TEL_ID = 21
TEL = f"tel_{TEL_ID:03d}"
# diffuse gammas, pointing at 70 deg altitude
SIMTEL_FILE = TEST_DATA_DIR / "simtel" / "gamma_200_800E3GeV_20_20deg_ATM52_102168.corsika.gz.NSBmed4.simtel.gz"
POINTING_ALT, POINTING_AZ = 70, 175.196


@pytest.fixture(scope="module")
def observation_dl1(observation_monitoring, tmp_path_factory):
    """DL1 file of the observation of tel 21, with its monitoring data"""
    run = OBSERVATION_RUNS[TEL_ID]
    output = tmp_path_factory.mktemp("dl1") / "observation.dl1.h5"
    run_process(run["path"], output, f"--max-events={run['n_events']}", *monitoring_source_args(observation_monitoring[TEL_ID]))
    return output


@pytest.fixture(scope="module")
def simulation_dl1(tmp_path_factory):
    output = tmp_path_factory.mktemp("dl1") / "gamma.dl1.h5"
    run_process(SIMTEL_FILE, output, config=MC_CONFIG_FILES["low"])
    return output


def test_observation_parameters(observation_dl1):
    data = load_dl1_sst1m(str(observation_dl1), tel=TEL)
    parameters = read_table(observation_dl1, f"/dl1/event/telescope/parameters/{TEL}")
    np.testing.assert_array_equal(data["event_id"], parameters["event_id"])

    intensity = np.asarray(parameters["camera_frame_hillas_intensity"])
    np.testing.assert_allclose(data["log_camera_frame_hillas_intensity"], np.log10(intensity))
    np.testing.assert_allclose(
        data["camera_frame_hillas_wl"],
        np.asarray(parameters["camera_frame_hillas_width"]) / np.asarray(parameters["camera_frame_hillas_length"]),
    )
    assert np.all(data["equivalent_focal_length"] == SubarrayDescription.from_hdf(observation_dl1).tel[TEL_ID].optics.equivalent_focal_length.to_value(u.m))

    # unix time of the triggers, date of the run (obs_id 202601201179)
    trigger_time = read_trigger_time(observation_dl1, TEL)
    np.testing.assert_array_equal(data["local_time"], trigger_time["trigger_time"].unix)
    assert set(data["date"]) == {20260120}
    # mean charge of the pedestal images of the single chunk of pedestal events
    image = read_table(observation_dl1, f"{DL1_PIXEL_STATISTICS_GROUP}/sky_pedestal_image/{TEL}")
    assert len(image) == 1
    np.testing.assert_allclose(data["meanQ"], np.nanmean(image["mean"][0]))
    assert 0 < data["meanQ"][0] < 10
    assert np.all((74 < data["true_alt_tel"]) & (data["true_alt_tel"] < 75))

    # the parameters are also in the pandas table
    df = load_dl1_sst1m(str(observation_dl1), tel=TEL, table="pandas")
    assert {"log_camera_frame_hillas_intensity", "local_time", "date", "meanQ"} <= set(df.columns)


def test_scaled_intensity(observation_dl1):
    """the log of the intensity is computed from the scaled intensity"""
    config = deepcopy(load_config(DATA_CONFIG_FILE, ismc=False))
    config["NsbCalibrator"]["intensity_correction"][TEL] = 2.0
    data = load_dl1_sst1m(str(observation_dl1), tel=TEL, config=config, scale_intensities=True)
    intensity = np.asarray(read_table(observation_dl1, f"/dl1/event/telescope/parameters/{TEL}")["camera_frame_hillas_intensity"])
    np.testing.assert_allclose(data["log_camera_frame_hillas_intensity"], np.log10(2 * intensity))


def test_simulation_parameters(simulation_dl1):
    tel = "tel_001"
    data = load_dl1_sst1m(str(simulation_dl1), tel=tel)
    parameters = read_table(simulation_dl1, f"/dl1/event/telescope/parameters/{tel}")
    np.testing.assert_array_equal(data["event_id"], parameters["event_id"])

    # pointing of the simulation
    np.testing.assert_allclose(data["true_alt_tel"], POINTING_ALT)
    np.testing.assert_allclose(data["true_az_tel"], POINTING_AZ, atol=1e-3)
    # true parameters of the showers
    shower = read_table(simulation_dl1, "/simulation/event/subarray/shower")
    true_energy = dict(zip(shower["event_id"], shower["true_energy"].quantity.to_value(u.TeV), strict=True))
    np.testing.assert_allclose(data["true_energy"], [true_energy[e] for e in data["event_id"]])
    np.testing.assert_allclose(data["log_true_energy"], np.log10(np.asarray(data["true_energy"])))

    # DISP (lstchain definition): vector from the center of gravity of the image to the true
    # position in the camera, its norm is its projection on the major axis
    image = np.isfinite(np.asarray(data["camera_frame_hillas_intensity"]))
    assert image.sum() > 50
    dx = np.asarray(data["true_camera_x"]) - np.asarray(data["camera_frame_hillas_x"])
    dy = np.asarray(data["true_camera_y"]) - np.asarray(data["camera_frame_hillas_y"])
    psi = np.deg2rad(np.asarray(data["camera_frame_hillas_psi"]))
    np.testing.assert_allclose(np.asarray(data["disp_dx"])[image], dx[image])
    np.testing.assert_allclose(np.asarray(data["disp_norm"])[image], np.abs(dx * np.cos(psi) + dy * np.sin(psi))[image])
    assert set(np.asarray(data["disp_sign"])[image]) <= {-1, 1}
    assert np.all(np.abs(np.asarray(data["miss"])[image]) < np.asarray(data["disp_norm"])[image] + 1e-9)
    # true impact distance: distance between the telescope and the core of the shower
    position = SubarrayDescription.from_hdf(simulation_dl1).positions[1]
    np.testing.assert_allclose(
        data["true_tel_impact_distance"],
        np.hypot(np.asarray(data["true_core_x"]) - position[0].to_value(u.m), np.asarray(data["true_core_y"]) - position[1].to_value(u.m)),
    )
