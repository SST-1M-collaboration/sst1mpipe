"""Machine learning reconstruction with the ctapipe reconstructors: training and application tools"""

import astropy.units as u
import numpy as np
import pytest
import tables
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, angular_separation
from astropy.time import Time
from astropy.table import Table, join
from ctapipe.containers import CoordinateFrameType
from ctapipe.core import run_tool
from ctapipe.io import read_table
from ctapipe.reco import Reconstructor

from conftest import OBSERVATION_RUNS, monitoring_source_args, run_process
from sst1mpipe.scripts.sst1mpipe_apply_models_tool import ApplyModels, add_telescope_pointing
from sst1mpipe.scripts.sst1mpipe_train_models_tool import (
    TrainDispReconstructor,
    TrainEnergyRegressor,
    TrainParticleClassifier,
)
from sst1mpipe.resources import (
    DATA_CONFIG_FILE,
    MC_CONFIG_FILES,
    STEREO_TRAIN_CONFIG_FILES,
    TEST_DATA_DIR,
    TRAIN_CONFIG_FILES,
)

SIMTEL_DIR = TEST_DATA_DIR / "simtel"
GAMMA_DIFFUSE = SIMTEL_DIR / "gamma_200_800E3GeV_20_20deg_ATM52_102168.corsika.gz.NSBmed4.simtel.gz"
GAMMA_POINT = SIMTEL_DIR / "gamma_point_200_800E3GeV_20_20deg_ATM52_300669.corsika.gz.NSBmed4.simtel.gz"
PROTON = SIMTEL_DIR / "proton_400_1300E3GeV_20_20deg_ATM52_204508.corsika.gz.NSBmed4.simtel.gz"


def test_add_telescope_pointing():
    """pointing of the target, in alt/az at the time of the events, without pointing of the telescope"""
    location = EarthLocation(lon=14.7825 * u.deg, lat=49.9124 * u.deg, height=511 * u.m)
    time = Time(["2026-01-20T22:00:00", "2026-01-20T23:00:00"])
    table = Table(dict(time=time, subarray_pointing_lon=[166.99, 166.99] * u.deg,
                       subarray_pointing_lat=[38.11, 38.11] * u.deg,
                       subarray_pointing_frame=[CoordinateFrameType.ICRS.value] * 2))
    assert add_telescope_pointing(table, location)
    expected = SkyCoord(ra=166.99 * u.deg, dec=38.11 * u.deg).transform_to(AltAz(obstime=time, location=location))
    np.testing.assert_allclose(table["telescope_pointing_altitude"].quantity.to_value(u.deg), expected.alt.deg)
    np.testing.assert_allclose(table["telescope_pointing_azimuth"].quantity.to_value(u.deg), expected.az.deg)
    # the pointing of the telescope is kept
    assert not add_telescope_pointing(table, location)
    # pointing of the target in alt/az (simulations)
    altaz = Table(dict(time=time, subarray_pointing_lon=[180, 180] * u.deg, subarray_pointing_lat=[70, 70] * u.deg,
                       subarray_pointing_frame=[CoordinateFrameType.ALTAZ.value] * 2))
    assert add_telescope_pointing(altaz, location)
    np.testing.assert_allclose(altaz["telescope_pointing_altitude"].quantity.to_value(u.deg), 70)


@pytest.fixture(scope="module")
def dl1_files(tmp_path_factory):
    """DL1 files of the test simulations, with the telescope frame parameters"""
    directory = tmp_path_factory.mktemp("ml_dl1")
    files = {}
    for name, path in (("gamma", GAMMA_DIFFUSE), ("gamma_point", GAMMA_POINT), ("proton", PROTON)):
        files[name] = directory / f"{name}.dl1.h5"
        run_process(path, files[name], config=MC_CONFIG_FILES["low"])
    return files


@pytest.fixture(scope="module")
def models(dl1_files, tmp_path_factory):
    """the three models of an NSB level, trained with the sst1mpipe tools"""
    directory = tmp_path_factory.mktemp("ml_models")
    node = directory / "nsb100"
    node.mkdir()
    tools = {}
    for name, tool, inputs in (
        ("energy_regressor", TrainEnergyRegressor(), ["-i", dl1_files["gamma"]]),
        ("particle_classifier", TrainParticleClassifier(),
         ["--signal", dl1_files["gamma"], "--background", dl1_files["proton"]]),
        ("disp_reconstructor", TrainDispReconstructor(), ["-i", dl1_files["gamma"]]),
    ):
        argv = [*map(str, inputs), "-o", str(node / f"{name}.pkl"), f"--config={TRAIN_CONFIG_FILES[name]}",
                "--CrossValidator.n_cross_validations=2"]
        assert run_tool(tool, argv=argv, raises=True) == 0
        tools[name] = tool
    return directory, tools


@pytest.mark.parametrize("name", ["energy_regressor", "particle_classifier", "disp_reconstructor"])
def test_train_tools_plots(models, name):
    """the models, their cross validation and their plots (feature importances, performance)"""
    directory, tools = models
    node = directory / "nsb100"
    assert (node / f"{name}.pkl").exists()
    plot_dir = node / f"{name}_plots"
    # one figure per quantity and per telescope
    quantities = {
        "energy_regressor": ["feature_importance", "resolution", "bias", "migration"],
        "particle_classifier": ["feature_importance", "gammaness", "roc", "roc_per_energy"],
        "disp_reconstructor": ["feature_importance_norm", "feature_importance_sign", "disp_migration",
                               "disp_resolution", "disp_bias", "disp_sign"],
    }[name]
    expected = {f"{name}_tel{tel_id}_{quantity}.png" for tel_id in (1, 2) for quantity in quantities}
    assert {p.name for p in tools[name].plots} == expected
    for path in tools[name].plots:
        assert path.parent == plot_dir and path.stat().st_size > 0
    model = {"energy_regressor": "regressor", "particle_classifier": "classifier", "disp_reconstructor": "models"}[name]
    assert len(getattr(tools[name], model)._models) == 2  # one model per telescope (two telescope types)


def test_apply_models(models, dl1_files, tmp_path):
    """the models (cos_zenith feature) are applied: ctapipe DL2, stereo direction"""
    directory, _ = models
    output = tmp_path / "gamma_point.dl2.h5"
    reconstructors = [arg for path in sorted((directory / "nsb100").glob("*.pkl")) for arg in ("-r", str(path))]
    argv = ["-i", str(dl1_files["gamma_point"]), *reconstructors, "-o", str(output), "--no-progress"]
    assert run_tool(ApplyModels(), argv=argv, raises=True) == 0

    with tables.open_file(output) as h5file:
        paths = {node._v_pathname for node in h5file.walk_nodes("/dl2/event/subarray", "Table")}
    assert {"/dl2/event/subarray/energy/RandomForestRegressor", "/dl2/event/subarray/geometry/disp",
            "/dl2/event/subarray/particle_type/RandomForestClassifier"} <= paths

    # stereo direction of the events (gammas from the pointing direction, 20 deg zenith)
    geometry = join(read_table(output, "/dl2/event/subarray/geometry/disp"),
                    read_table(output, "/simulation/event/subarray/shower"), keys=["obs_id", "event_id"])
    stereo = geometry["disp_is_valid"] & (np.sum(geometry["disp_telescopes"], axis=1) == 2)
    assert np.count_nonzero(stereo) > 50
    separation = angular_separation(
        geometry["disp_az"].quantity, geometry["disp_alt"].quantity,
        geometry["true_az"].quantity, geometry["true_alt"].quantity,
    )[stereo]
    assert np.median(separation) < 1 * u.deg


@pytest.mark.parametrize("with_monitoring", [True, False])
def test_apply_models_to_observations(models, observation_monitoring, tmp_path, with_monitoring):
    """
    the models of the simulated telescope 1 are applied to the observations of telescope 21; without
    monitoring data (no pointing of the telescope), with the pointing of the target (warning)
    """
    run = OBSERVATION_RUNS[21]
    dl1 = tmp_path / "observation.dl1.h5"
    monitoring = monitoring_source_args(observation_monitoring[21]) if with_monitoring else []
    run_process(run["path"], dl1, f"--max-events={run['n_events']}", *monitoring, config=DATA_CONFIG_FILE)
    directory, _ = models
    output = tmp_path / "observation.dl2.h5"
    reconstructors = [arg for path in sorted((directory / "nsb100").glob("*.pkl")) for arg in ("-r", str(path))]
    argv = ["-i", str(dl1), *reconstructors, "-o", str(output), "--no-progress"]
    tool = ApplyModels()
    assert run_tool(tool, argv=argv, raises=True) == 0
    assert getattr(tool, "_warned_pointing", False) is not with_monitoring
    energy = read_table(output, "/dl2/event/telescope/energy/RandomForestRegressor/tel_021")
    assert len(energy) > 0 and np.any(energy["RandomForestRegressor_tel_is_valid"])
    geometry = read_table(output, "/dl2/event/subarray/geometry/disp")
    assert np.any(geometry["disp_is_valid"])


def test_apply_models_one_telescope(models, dl1_files, tmp_path):
    """telescopes=[1]: the DL2 of the subarray is the mono reconstruction of telescope 1"""
    directory, _ = models
    output = tmp_path / "gamma_point_tel1.dl2.h5"
    reconstructors = [arg for path in sorted((directory / "nsb100").glob("*.pkl")) for arg in ("-r", str(path))]
    argv = ["-i", str(dl1_files["gamma_point"]), *reconstructors, "-o", str(output), "--no-progress",
            "--ApplyModels.telescopes=[1]"]
    assert run_tool(ApplyModels(), argv=argv, raises=True) == 0
    stereo = read_table(output, "/dl2/event/subarray/energy/RandomForestRegressor")
    # one row per triggered event, telescope 2 never valid
    assert len(stereo) == len(read_table(output, "/dl1/event/subarray/trigger"))
    tel2 = read_table(output, "/dl2/event/telescope/energy/RandomForestRegressor/tel_002")
    assert not np.any(tel2["RandomForestRegressor_tel_is_valid"])
    mono = read_table(output, "/dl2/event/telescope/energy/RandomForestRegressor/tel_001")
    events = join(stereo, mono, keys=["obs_id", "event_id"])
    valid = events["RandomForestRegressor_is_valid"]
    assert np.count_nonzero(valid) > 0
    np.testing.assert_allclose(events["RandomForestRegressor_energy"][valid],
                               events["RandomForestRegressor_tel_energy"][valid], rtol=1e-6)
    assert np.all(np.asarray(stereo["RandomForestRegressor_telescopes"])[:, 1] == 0)



def test_stereo_models(dl1_files, tmp_path):
    """
    stereo models: the stereo parameters of the HillasReconstructor (impact distance, h_max) as features,
    trained on the events reconstructed by both telescopes; predictions only for these events
    """
    directory = tmp_path / "stereo"
    directory.mkdir()
    for name, tool, inputs in (
        ("energy_regressor", TrainEnergyRegressor(), ["-i", dl1_files["gamma"]]),
        ("particle_classifier", TrainParticleClassifier(),
         ["--signal", dl1_files["gamma"], "--background", dl1_files["proton"]]),
        ("disp_reconstructor", TrainDispReconstructor(), ["-i", dl1_files["gamma"]]),
    ):
        argv = [*map(str, inputs), "-o", str(directory / f"{name}.pkl"),
                f"--config={STEREO_TRAIN_CONFIG_FILES[name]}", "--CrossValidator.n_cross_validations=2"]
        assert run_tool(tool, argv=argv, raises=True) == 0
        model = Reconstructor.read(directory / f"{name}.pkl")
        assert {"HillasReconstructor_tel_impact_distance", "HillasReconstructor_h_max"} <= set(model.features)

    output = tmp_path / "gamma_point.dl2.h5"
    reconstructors = [arg for path in sorted(directory.glob("*.pkl")) for arg in ("-r", str(path))]
    argv = ["-i", str(dl1_files["gamma_point"]), *reconstructors, "-o", str(output), "--no-progress"]
    assert run_tool(ApplyModels(), argv=argv, raises=True) == 0
    energy = read_table(output, "/dl2/event/subarray/energy/RandomForestRegressor")
    geometry = read_table(output, "/dl2/event/subarray/geometry/HillasReconstructor")
    events = join(energy, geometry, keys=["obs_id", "event_id"])
    valid = np.asarray(events["RandomForestRegressor_is_valid"], dtype=bool)
    assert np.any(valid)
    # only the events with the stereo parameters
    assert np.all(np.isfinite(events["HillasReconstructor_h_max"][valid]))
