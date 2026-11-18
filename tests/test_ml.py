"""Machine learning reconstruction with the ctapipe reconstructors (sst1mpipe.ml): MARS stereo combination, tools"""

import astropy.units as u
import numpy as np
import pytest
import tables
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, angular_separation
from astropy.time import Time
from astropy.table import Table, join
from ctapipe.containers import (
    CoordinateFrameType,
    ArrayEventContainer,
    DispContainer,
    ImageParametersContainer,
    ReconstructedGeometryContainer,
)
from ctapipe.core import run_tool
from ctapipe.io import read_table
from ctapipe.reco.preprocessing import telescope_to_horizontal
from ctapipe.reco.reconstructor import ReconstructionProperty
from ctapipe.reco.stereo_combination import StereoCombiner, StereoMeanCombiner

from conftest import run_process
from sst1mpipe.ml import StereoMARSCombiner, mars_disp_signs
from sst1mpipe.scripts.sst1mpipe_apply_models_tool import ApplyModels, add_telescope_pointing
from sst1mpipe.scripts.sst1mpipe_train_models_tool import (
    TrainDispReconstructor,
    TrainEnergyRegressor,
    TrainParticleClassifier,
)
from sst1mpipe.resources import MC_CONFIG_FILES, TEST_DATA_DIR, TRAIN_CONFIG_FILES

SIMTEL_DIR = TEST_DATA_DIR / "simtel"
GAMMA_DIFFUSE = SIMTEL_DIR / "gamma_200_800E3GeV_20_20deg_ATM52_102168.corsika.gz.NSBmed4.simtel.gz"
GAMMA_POINT = SIMTEL_DIR / "gamma_point_200_800E3GeV_20_20deg_ATM52_300669.corsika.gz.NSBmed4.simtel.gz"
PROTON = SIMTEL_DIR / "proton_400_1300E3GeV_20_20deg_ATM52_204508.corsika.gz.NSBmed4.simtel.gz"


def test_stereo_mars_combiner_is_available():
    """the combiner can be chosen by name in the configuration of the reconstructors"""
    assert "StereoMARSCombiner" in StereoCombiner.non_abstract_subclasses()


def test_mars_disp_signs():
    """each telescope takes the sign of its disp whose position is the closest to the crossing of the axes"""
    # major axes crossing at (0, 0): tel 1 centroid (1, 0) along lon, tel 2 centroid (0, 2) along lat
    hillas = [dict(fov_lon=1 * u.deg, fov_lat=0 * u.deg, psi=0 * u.deg),
              dict(fov_lon=0 * u.deg, fov_lat=2 * u.deg, psi=90 * u.deg)]
    for disps in ([0.9, 1.5], [-0.9, 1.5], [-0.9, -1.5]):
        np.testing.assert_array_equal(mars_disp_signs(hillas, disps * u.deg), [-1, -1])
    # crossing on the other side of the centroids
    hillas[0]["fov_lon"], hillas[1]["fov_lat"] = -1 * u.deg, -2 * u.deg
    np.testing.assert_array_equal(mars_disp_signs(hillas, [-0.9, -1.5] * u.deg), [1, 1])
    # parallel axes: no crossing, the predicted signs
    parallel = [dict(fov_lon=0, fov_lat=0, psi=30), dict(fov_lon=0, fov_lat=1, psi=30)]
    np.testing.assert_array_equal(mars_disp_signs(parallel, [0.5, -0.5]), [1, -1])

    with pytest.raises(ValueError, match="two telescopes"):
        mars_disp_signs(hillas[:1], [0.9])
    with pytest.raises(NotImplementedError, match="more than two telescopes"):
        mars_disp_signs([*hillas, hillas[0]], [0.9, 1.5, 0.3])


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
    table["disp_tel_alt"], table["disp_tel_az"] = direction(table, table["disp_tel_parameter"].quantity)
    return table


def direction(table, disp):
    psi = table["hillas_psi"].quantity
    alt, az = telescope_to_horizontal(
        lon=table["hillas_fov_lon"].quantity + disp * np.cos(psi),
        lat=table["hillas_fov_lat"].quantity + disp * np.sin(psi),
        pointing_alt=table["telescope_pointing_altitude"].quantity,
        pointing_az=table["telescope_pointing_azimuth"].quantity,
    )
    return alt.to(u.deg), az.to(u.deg)


def test_stereo_mars_combiner_table():
    """the StereoMeanCombiner of the directions with the MARS disp signs (events with two telescopes)"""
    table = mono_table()
    stereo = StereoMARSCombiner(prefix="disp", property=ReconstructionProperty.GEOMETRY,
                                weights="intensity").predict_table(table)

    corrected = table.copy()
    for event_id in np.unique(table["event_id"]):
        rows = np.flatnonzero((table["event_id"] == event_id) & table["disp_tel_is_valid"])
        if len(rows) == 2:
            hillas = [{name: table[f"hillas_{name}"].quantity[row] for name in ("fov_lon", "fov_lat", "psi")}
                      for row in rows]
            disp = table["disp_tel_parameter"].quantity[rows]
            corrected["disp_tel_parameter"][rows] = (mars_disp_signs(hillas, disp) * np.abs(disp)).to_value(u.deg)
    corrected["disp_tel_alt"], corrected["disp_tel_az"] = direction(corrected, corrected["disp_tel_parameter"].quantity)
    expected = StereoMeanCombiner(prefix="disp", property=ReconstructionProperty.GEOMETRY,
                                  weights="intensity").predict_table(corrected)

    np.testing.assert_array_equal(stereo["disp_is_valid"], expected["disp_is_valid"])
    valid = stereo["disp_is_valid"]
    separation = angular_separation(stereo["disp_az"].quantity, stereo["disp_alt"].quantity,
                                    expected["disp_az"].quantity, expected["disp_alt"].quantity)[valid]
    assert np.all(separation < 1e-6 * u.deg)
    # at least one sign corrected
    assert np.any(corrected["disp_tel_parameter"] != table["disp_tel_parameter"])


def test_stereo_mars_combiner_more_than_two_telescopes():
    table = mono_table(n_events=2, n_tel=3)
    table["disp_tel_is_valid"] = True
    combiner = StereoMARSCombiner(prefix="disp", property=ReconstructionProperty.GEOMETRY)
    with pytest.raises(NotImplementedError):
        combiner.predict_table(table)


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
        event.dl2.tel[tel_id].geometry["disp"] = ReconstructedGeometryContainer(
            alt=row["disp_tel_alt"] * u.deg, az=row["disp_tel_az"] * u.deg, is_valid=True)
    combiner(event)
    stereo = event.dl2.stereo.geometry["disp"]
    assert stereo.is_valid and stereo.telescopes == [1, 2]
    separation = angular_separation(stereo.az, stereo.alt, expected["disp_az"] * u.deg, expected["disp_alt"] * u.deg)
    assert separation < 1e-6 * u.deg
    # the directions of the telescopes keep their predicted signs
    for row in table:
        assert event.dl2.tel[row["tel_id"]].geometry["disp"].alt == row["disp_tel_alt"] * u.deg


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
    assert sorted(p.name for p in tools[name].plots) == [f"{name}_feature_importance.png", f"{name}_performance.png"]
    for path in tools[name].plots:
        assert path.parent == plot_dir and path.stat().st_size > 0
    model = {"energy_regressor": "regressor", "particle_classifier": "classifier", "disp_reconstructor": "models"}[name]
    assert len(getattr(tools[name], model)._models) == 2  # one model per telescope (two telescope types)


def test_apply_models(models, dl1_files, tmp_path):
    """the models (cos_zenith feature) are applied: ctapipe DL2, MARS stereo direction"""
    directory, _ = models
    output = tmp_path / "gamma_point.dl2.h5"
    reconstructors = [arg for path in sorted((directory / "nsb100").glob("*.pkl")) for arg in ("-r", str(path))]
    argv = ["-i", str(dl1_files["gamma_point"]), *reconstructors, "-o", str(output), "--no-progress"]
    assert run_tool(ApplyModels(), argv=argv, raises=True) == 0

    with tables.open_file(output) as h5file:
        paths = {node._v_pathname for node in h5file.walk_nodes("/dl2/event/subarray", "Table")}
    assert {"/dl2/event/subarray/energy/RandomForestRegressor", "/dl2/event/subarray/geometry/disp",
            "/dl2/event/subarray/particle_type/RandomForestClassifier"} <= paths

    # MARS direction of the stereo events (gammas from the pointing direction, 20 deg zenith)
    geometry = join(read_table(output, "/dl2/event/subarray/geometry/disp"),
                    read_table(output, "/simulation/event/subarray/shower"), keys=["obs_id", "event_id"])
    stereo = geometry["disp_is_valid"] & (np.sum(geometry["disp_telescopes"], axis=1) == 2)
    assert np.count_nonzero(stereo) > 50
    separation = angular_separation(
        geometry["disp_az"].quantity, geometry["disp_alt"].quantity,
        geometry["true_az"].quantity, geometry["true_alt"].quantity,
    )[stereo]
    assert np.median(separation) < 1 * u.deg
