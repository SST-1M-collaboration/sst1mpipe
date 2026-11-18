"""Machine learning reconstruction with the ctapipe reconstructors (sst1mpipe.ml): MARS stereo combination, tools"""
import itertools

import astropy.units as u
import numpy as np
import pytest
import tables
from astropy.coordinates import angular_separation
from astropy.table import Table, join
from ctapipe.containers import (
    ArrayEventContainer,
    DispContainer,
    ImageParametersContainer,
    ReconstructedGeometryContainer,
)
from ctapipe.core import run_tool
from ctapipe.io import read_table
from ctapipe.reco.preprocessing import telescope_to_horizontal
from ctapipe.reco.reconstructor import ReconstructionProperty
from ctapipe.reco.stereo_combination import StereoCombiner

from conftest import OBSERVATION_RUNS, monitoring_source_args, run_process
from sst1mpipe.ml import StereoMARSCombiner, mars_directions
from sst1mpipe.scripts.sst1mpipe_apply_models_tool import ApplyModels, closest_model_node, find_model_nodes
from sst1mpipe.scripts.sst1mpipe_train_models_tool import (
    TrainDispReconstructor,
    TrainEnergyRegressor,
    TrainParticleClassifier,
)
from sst1mpipe.resources import DATA_CONFIG_FILE, MC_CONFIG_FILES, TEST_DATA_DIR, TRAIN_CONFIG_FILES

SIMTEL_DIR = TEST_DATA_DIR / "simtel"
GAMMA_DIFFUSE = SIMTEL_DIR / "gamma_200_800E3GeV_20_20deg_ATM52_102168.corsika.gz.NSBmed4.simtel.gz"
GAMMA_POINT = SIMTEL_DIR / "gamma_point_200_800E3GeV_20_20deg_ATM52_300669.corsika.gz.NSBmed4.simtel.gz"
PROTON = SIMTEL_DIR / "proton_400_1300E3GeV_20_20deg_ATM52_204508.corsika.gz.NSBmed4.simtel.gz"


def unit_vector(alt, az):
    return np.array([np.cos(alt) * np.cos(az), np.cos(alt) * np.sin(az), np.sin(alt)])


def test_stereo_mars_combiner_is_available():
    """the combiner can be chosen by name in the configuration of the reconstructors"""
    assert "StereoMARSCombiner" in StereoCombiner.non_abstract_subclasses()


def test_mars_directions():
    """the combination of signs whose positions are the closest is chosen, then weighted"""
    alt = np.deg2rad(70)
    # telescope 1: predicted position at az 0.0, other at az 1.0; telescope 2: predicted 2.0, other 1.02
    candidates = np.array([
        [unit_vector(alt, np.deg2rad(0.0)), unit_vector(alt, np.deg2rad(1.0))],
        [unit_vector(alt, np.deg2rad(2.0)), unit_vector(alt, np.deg2rad(1.02))],
        # array event 1: a single telescope keeps its predicted position
        [unit_vector(alt, np.deg2rad(5.0)), unit_vector(alt, np.deg2rad(6.0))],
    ])
    weights = np.array([1.0, 3.0, 1.0])
    alt_reco, az_reco, spread, distance = mars_directions(candidates, weights, np.array([0, 0, 1]), 3)
    mean = candidates[0, 1] + 3 * candidates[1, 1]
    np.testing.assert_allclose(alt_reco[:2], np.arcsin([mean[2] / np.linalg.norm(mean), candidates[2, 0, 2]]))
    np.testing.assert_allclose(az_reco[:2], np.deg2rad([(1.0 + 3 * 1.02) / 4, 5.0]), atol=1e-4)
    separation = np.arccos(np.dot(candidates[0, 1], candidates[1, 1]))
    np.testing.assert_allclose(distance[:2], [separation, 0], atol=1e-12)
    assert spread[1] < 1e-6
    # array event 2 without telescope event
    assert np.isnan(alt_reco[2]) and np.isnan(distance[2])


def mono_table(n_events=40, n_tel=2, seed=0):
    """telescope events with the columns of the disp reconstructor (pointing 70 deg altitude)"""
    rng = np.random.default_rng(seed)
    rows = []
    for event_id in range(n_events):
        for tel_id in range(1, n_tel + 1):
            rows.append(dict(
                obs_id=1, event_id=event_id, tel_id=tel_id,
                hillas_fov_lon=rng.uniform(-1, 1), hillas_fov_lat=rng.uniform(-1, 1),
                hillas_psi=rng.uniform(-90, 90), hillas_intensity=rng.uniform(50, 5000),
                disp_tel_parameter=rng.uniform(-1.5, 1.5), disp_tel_is_valid=rng.uniform() > 0.1,
                telescope_pointing_altitude=70.0, telescope_pointing_azimuth=10.0,
            ))
    table = Table(rows=rows)
    for name in ("hillas_fov_lon", "hillas_fov_lat", "hillas_psi", "disp_tel_parameter",
                 "telescope_pointing_altitude", "telescope_pointing_azimuth"):
        table[name].unit = u.deg
    return table


def brute_force(table, weights="intensity"):
    """stereo directions by explicit loops: (alt, az) in deg per event_id, NaN without valid telescope"""
    result = {}
    for event_id in np.unique(table["event_id"]):
        rows = table[(table["event_id"] == event_id) & table["disp_tel_is_valid"]]
        if len(rows) == 0:
            result[event_id] = (np.nan, np.nan)
            continue
        positions = []
        for row in rows:
            candidates = []
            for sign in (1, -1):
                disp = sign * row["disp_tel_parameter"]
                psi = np.deg2rad(row["hillas_psi"])
                alt, az = telescope_to_horizontal(
                    lon=(row["hillas_fov_lon"] + disp * np.cos(psi)) * u.deg,
                    lat=(row["hillas_fov_lat"] + disp * np.sin(psi)) * u.deg,
                    pointing_alt=row["telescope_pointing_altitude"] * u.deg,
                    pointing_az=row["telescope_pointing_azimuth"] * u.deg,
                )
                candidates.append((alt, az))
            positions.append(candidates)
        # MARS sign choice for two telescopes only, the predicted signs otherwise
        best, best_distance = [p[0] for p in positions], np.inf
        for signs in itertools.product((0, 1), repeat=len(rows)) if len(rows) == 2 else []:
            chosen = [positions[i][s] for i, s in enumerate(signs)]
            distance = sum(
                angular_separation(a[1], a[0], b[1], b[0]).to_value(u.rad)
                for a, b in itertools.combinations(chosen, 2)
            )
            if distance < best_distance:
                best, best_distance = chosen, distance
        w = np.asarray(rows["hillas_intensity"]) if weights == "intensity" else np.ones(len(rows))
        vector = sum(wi * unit_vector(alt.to_value(u.rad), az.to_value(u.rad)) for wi, (alt, az) in zip(w, best, strict=True))
        result[event_id] = (np.rad2deg(np.arcsin(vector[2] / np.linalg.norm(vector))),
                            np.rad2deg(np.arctan2(vector[1], vector[0])) % 360)
    return result


@pytest.mark.parametrize("n_tel", [2, 3])
def test_stereo_mars_combiner_table(n_tel):
    """stereo directions of a table, compared with explicit loops over the events and sign combinations"""
    table = mono_table(n_tel=n_tel)
    combiner = StereoMARSCombiner(prefix="disp", property=ReconstructionProperty.GEOMETRY, weights="intensity")
    stereo = combiner.predict_table(table)
    expected = brute_force(table)
    assert len(stereo) == len(expected)
    for row in stereo:
        alt, az = expected[row["event_id"]]
        if np.isnan(alt):
            assert not row["disp_is_valid"]
            continue
        assert row["disp_is_valid"]
        separation = angular_separation(row["disp_az"] * u.deg, row["disp_alt"] * u.deg, az * u.deg, alt * u.deg)
        assert separation < 1e-6 * u.deg
        assert np.isfinite(row["disp_goodness_of_fit"])
    np.testing.assert_array_equal(
        [len(t) for t in stereo["disp_telescopes"]],
        [np.count_nonzero(table["disp_tel_is_valid"][table["event_id"] == e]) for e in stereo["event_id"]],
    )


def test_stereo_mars_combiner_event():
    """the event loop version gives the direction of the table version"""
    table = mono_table(n_events=1)
    table["disp_tel_is_valid"] = True
    combiner = StereoMARSCombiner(prefix="disp", property=ReconstructionProperty.GEOMETRY, weights="none")
    expected = combiner.predict_table(table)[0]

    event = ArrayEventContainer()
    for row in table:
        tel_id = row["tel_id"]
        event.dl1.tel[tel_id].parameters = ImageParametersContainer()
        hillas = event.dl1.tel[tel_id].parameters.hillas
        hillas.fov_lon, hillas.fov_lat = row["hillas_fov_lon"] * u.deg, row["hillas_fov_lat"] * u.deg
        hillas.psi, hillas.intensity = row["hillas_psi"] * u.deg, row["hillas_intensity"]
        event.monitoring.tel[tel_id].pointing.altitude = row["telescope_pointing_altitude"] * u.deg
        event.monitoring.tel[tel_id].pointing.azimuth = row["telescope_pointing_azimuth"] * u.deg
        event.dl2.tel[tel_id].disp["disp"] = DispContainer(parameter=row["disp_tel_parameter"] * u.deg)
        event.dl2.tel[tel_id].geometry["disp"] = ReconstructedGeometryContainer(is_valid=True)
    combiner(event)
    stereo = event.dl2.stereo.geometry["disp"]
    assert stereo.is_valid and stereo.telescopes == [1, 2]
    separation = angular_separation(stereo.az, stereo.alt, expected["disp_az"] * u.deg, expected["disp_alt"] * u.deg)
    assert separation < 1e-6 * u.deg
    assert np.isclose(stereo.goodness_of_fit, expected["disp_goodness_of_fit"])


def test_closest_model_node(tmp_path):
    for name in ("ze20_az0_nsb100", "ze20_az180_nsb100", "ze40_az0_nsb100", "ze20_az0_nsb200", "README"):
        (tmp_path / name).mkdir()
    nodes = find_model_nodes(tmp_path)
    assert len(nodes) == 4
    assert closest_model_node(nodes, zenith=22).name == "ze20_az0_nsb100"
    assert closest_model_node(nodes, zenith=35).name == "ze40_az0_nsb100"
    # azimuth through 360: 350 is closer to 0 than to 180
    assert closest_model_node(nodes, zenith=20, azimuth=350).name in ("ze20_az0_nsb100", "ze20_az0_nsb200")
    assert closest_model_node(nodes, zenith=20, azimuth=170).name == "ze20_az180_nsb100"
    assert closest_model_node(nodes, zenith=20, nsb_rate=180).name == "ze20_az0_nsb200"
    with pytest.raises(ValueError, match="No model node"):
        closest_model_node({}, zenith=20)


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
    """the three models of the node ze20_az0_nsb100 of a production of models, trained with the sst1mpipe tools"""
    directory = tmp_path_factory.mktemp("ml_models")
    node = directory / "ze20_az0_nsb100"
    node.mkdir()
    (directory / "ze50_az0_nsb100").mkdir()
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
    node = directory / "ze20_az0_nsb100"
    assert (node / f"{name}.pkl").exists()
    plot_dir = node / f"{name}_plots"
    assert sorted(p.name for p in tools[name].plots) == [f"{name}_feature_importance.png", f"{name}_performance.png"]
    for path in tools[name].plots:
        assert path.parent == plot_dir and path.stat().st_size > 0
    assert len(tools[name].model._models) == 2  # one model per telescope (two telescope types)


def test_apply_models(models, dl1_files, tmp_path):
    """the models of the closest node are applied: ctapipe DL2, MARS stereo direction"""
    directory, _ = models
    output = tmp_path / "gamma_point.dl2.h5"
    tool = ApplyModels()
    argv = ["-i", str(dl1_files["gamma_point"]), "--models-directory", str(directory), "-o", str(output), "--no-progress"]
    assert run_tool(tool, argv=argv, raises=True) == 0
    assert {p.parent.name for p in tool.reconstructor_paths} == {"ze20_az0_nsb100"}
    assert len(tool.reconstructor_paths) == 3

    with tables.open_file(output) as h5file:
        paths = {node._v_pathname for node in h5file.walk_nodes("/dl2/event/subarray", "Table")}
    assert {"/dl2/event/subarray/energy/RandomForestRegressor", "/dl2/event/subarray/geometry/disp",
            "/dl2/event/subarray/particle_type/RandomForestClassifier"} <= paths

    # MARS direction of the stereo events (gammas from the pointing direction, 20 deg zenith)
    geometry = join(read_table(output, "/dl2/event/subarray/geometry/disp"),
                    read_table(output, "/simulation/event/subarray/shower"), keys=["obs_id", "event_id"])
    stereo = geometry["disp_is_valid"] & (np.sum(geometry["disp_telescopes"], axis=1) == 2)
    assert np.count_nonzero(stereo) > 50
    assert np.all(np.isfinite(geometry["disp_goodness_of_fit"][stereo]))
    separation = angular_separation(
        geometry["disp_az"].quantity, geometry["disp_alt"].quantity,
        geometry["true_az"].quantity, geometry["true_alt"].quantity,
    )[stereo]
    assert np.median(separation) < 1 * u.deg


def test_apply_models_to_observations(models, observation_monitoring, tmp_path):
    """the models of the simulated telescope 1 are applied to the observations of telescope 21"""
    run = OBSERVATION_RUNS[21]
    dl1 = tmp_path / "observation.dl1.h5"
    run_process(run["path"], dl1, f"--max-events={run['n_events']}", *monitoring_source_args(observation_monitoring[21]),
                config=DATA_CONFIG_FILE)
    directory, _ = models
    output = tmp_path / "observation.dl2.h5"
    tool = ApplyModels()
    argv = ["-i", str(dl1), "--models-directory", str(directory), "-o", str(output), "--no-progress"]
    assert run_tool(tool, argv=argv, raises=True) == 0
    energy = read_table(output, "/dl2/event/telescope/energy/RandomForestRegressor/tel_021")
    assert len(energy) > 0 and np.any(energy["RandomForestRegressor_tel_is_valid"])
    geometry = read_table(output, "/dl2/event/subarray/geometry/disp")
    assert np.any(geometry["disp_is_valid"])

