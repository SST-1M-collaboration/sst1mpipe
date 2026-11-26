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



@pytest.fixture(scope="module")
def irfs(models, dl1_files, tmp_path_factory):
    """
    IRFs (cuts of the Crab paper: few events) of the DL2 of the models with the sst1mpipe IRF tools: (IRF
    directory: cuts.fits, irf_point_like.fits.gz, benchmarks_point_like.fits.gz; DL2 directory). The spectra
    of the ctapipe IRF tools (changed by the sst1mpipe tools) are restored after the tests.
    """
    from astropy.table import QTable
    from ctapipe.irf import spectra
    from pyirf.spectral import CRAB_MAGIC_JHEAP2015, DAMPE_P_He_SPECTRUM

    from sst1mpipe.resources import CONFIG_DIR
    from sst1mpipe.scripts.sst1mpipe_irf_tools import SST1MEventSelectionOptimizer, SST1MIrfTool

    saved_spectra = dict(spectra.SPECTRA)

    tmp_path = tmp_path_factory.mktemp("irfs")
    directory, _ = models
    reconstructors = [arg for path in sorted((directory / "nsb100").glob("*.pkl")) for arg in ("-r", str(path))]
    dl2 = tmp_path / "dl2"
    dl2.mkdir()
    for name, particle in (("gamma_point", "gamma_point"), ("proton_test", "proton")):
        argv = ["-i", str(dl1_files[particle]), *reconstructors, "-o", str(dl2 / f"{name}.dl2.h5"), "--no-progress"]
        assert run_tool(ApplyModels(), argv=argv, raises=True) == 0

    irf = tmp_path / "irf"
    irf.mkdir()
    configs = [f"--config={CONFIG_DIR / 'sst1mpipe_optimize_cuts.yaml'}",
               f"--config={CONFIG_DIR / 'sst1mpipe_irf_paper_cuts.yaml'}"]
    argv = [*configs, "-g", str(dl2 / "gamma_point.dl2.h5"), "-p", str(dl2 / "proton_test.dl2.h5"),
            "-o", str(irf / "cuts.fits")]
    assert run_tool(SST1MEventSelectionOptimizer(), argv=argv, raises=True) == 0
    # spectra of the paper: Crab of MAGIC, protons + helium
    assert spectra.SPECTRA[spectra.Spectra.CRAB_HEGRA] is CRAB_MAGIC_JHEAP2015
    assert spectra.SPECTRA[spectra.Spectra.IRFDOC_PROTON_SPECTRUM] is DAMPE_P_He_SPECTRUM
    # one theta cut for all the energies (GlobalThetaPercentileCuts), gammaness cut per energy bin
    rad_max = QTable.read(irf / "cuts.fits", hdu="RAD_MAX")["cut"]
    assert len(rad_max) > 1 and np.all(rad_max == rad_max[0])
    # FixedThetaCuts: the same gammaness cuts, a fixed theta cut (0.2 deg by default)
    import astropy.units as u

    for theta_cut, options in ((0.2, []), (0.25, ["--FixedThetaCuts.theta_cut=0.25 deg"])):
        output = irf / f"cuts_fixed_{theta_cut}.fits"
        argv = [*configs, "--EventSelectionOptimizer.optimization_algorithm=FixedThetaCuts", *options,
                "-g", str(dl2 / "gamma_point.dl2.h5"), "-p", str(dl2 / "proton_test.dl2.h5"), "-o", str(output)]
        assert run_tool(SST1MEventSelectionOptimizer(), argv=argv, raises=True) == 0
        fixed = QTable.read(output, hdu="RAD_MAX")["cut"]
        assert len(fixed) == len(rad_max) and np.allclose(fixed.to_value(u.deg), theta_cut)
        np.testing.assert_array_equal(QTable.read(output, hdu="GH_CUTS")["cut"], QTable.read(irf / "cuts.fits", hdu="GH_CUTS")["cut"])
    argv = [f"--config={CONFIG_DIR / 'sst1mpipe_compute_irf.yaml'}", "--cuts", str(irf / "cuts.fits"),
            "--gamma-file", str(dl2 / "gamma_point.dl2.h5"), "--proton-file", str(dl2 / "proton_test.dl2.h5"),
            "--do-background", "--spatial-selection-applied", "--DefaultFoVOffsetBins.fov_offset_n_bins=1",
            "--output", str(irf / "irf_point_like.fits.gz"), "--benchmark-output", str(irf / "benchmarks_point_like.fits.gz")]
    assert run_tool(SST1MIrfTool(), argv=argv, raises=True) == 0
    # background of the sensitivity within 1.2 deg (Sensitivity2dMaker), the IRFs within 5 deg (DefaultFoVOffsetBins)
    sensitivity = QTable.read(irf / "benchmarks_point_like.fits.gz", hdu="SENSITIVITY")
    assert u.Quantity(sensitivity["THETA_LO"]).ravel()[0] == 0 * u.deg
    assert u.Quantity(sensitivity["THETA_HI"]).ravel()[-1] == 1.2 * u.deg
    background = QTable.read(irf / "irf_point_like.fits.gz", hdu="BACKGROUND")
    assert u.Quantity(background["THETA_HI"]).ravel()[-1] == 5 * u.deg
    yield irf, dl2
    spectra.SPECTRA.update(saved_spectra)


def test_plot_irfs(irfs, tmp_path):
    """Figures of the IRFs, with references"""
    from sst1mpipe.scripts.sst1mpipe_plot_irfs_tool import main as plot_irfs

    irf, dl2 = irfs
    plot_irfs(["--irf", f"stereo={irf}", "--dl2", f"stereo={dl2}", "-o", str(tmp_path / "plots")])
    names = {path.stem for path in (tmp_path / "plots").glob("*.png")}
    assert names == {"effective_area", "effective_area_protons", "energy_bias", "energy_resolution", "angular_resolution",
                     "theta_cut", "gammaness_cut", "background", "sensitivity", "sensitivity_crab", "energy_migration",
                     "gammaness", "roc"}
    # comparison: one series per input
    plot_irfs(["--irf", f"a={irf}", "--irf", f"b={irf}", "-o", str(tmp_path / "comparison")])
    assert {"energy_migration_a", "energy_migration_b", "sensitivity"} <= {
        path.stem for path in (tmp_path / "comparison").glob("*.png")}
    # with the performance of the paper (IRF directories of the same format, some quantities only)
    plot_irfs(["--irf", f"stereo={irf}", "--reference", "paper 30=ref:crab_paper_2025/stereo/ze30",
               "--reference", "paper 20=ref:crab_paper_2025/stereo/ze20", "-o", str(tmp_path / "reference")])
    assert {path.stem for path in (tmp_path / "reference").glob("*.png")} == names - {"gammaness", "roc"}
    # the paper alone: its quantities only
    plot_irfs(["--reference", "tel1=ref:crab_paper_2025/tel1/ze30", "--reference", "tel2=ref:crab_paper_2025/tel2/ze20",
               "--reference", "ASTRI=ref:astri_50h", "-o", str(tmp_path / "paper")])
    assert {path.stem for path in (tmp_path / "paper").glob("*.png")} == {
        "effective_area", "energy_bias", "energy_resolution", "angular_resolution", "sensitivity", "sensitivity_crab"}


def test_crab_paper_reference():
    """Performance of the SST-1M Crab paper in IRF directories: sensible values"""
    import astropy.units as u

    from sst1mpipe.performance.irf_plots import REFERENCE_DIRECTORY, read_irf_directory

    stereo = read_irf_directory(REFERENCE_DIRECTORY / "crab_paper_2025" / "stereo" / "ze30")
    assert set(stereo) == {"energy", "angular", "sensitivity"}
    # best sensitivity of the stereo: ~2e-12 TeV cm-2 s-1 at ~10 TeV (Fig. 11)
    sensitivity = stereo["sensitivity"]["ENERGY_FLUX_SENSITIVITY"].to_value(u.TeV / u.cm**2 / u.s)
    assert 1.5e-12 < np.nanmin(sensitivity) < 2.5e-12
    aeff = read_irf_directory(REFERENCE_DIRECTORY / "crab_paper_2025" / "tel2" / "ze20")["aeff"]["EFFAREA"].to_value(u.m**2)
    assert 5e5 < np.nanmax(aeff) < 1e6


@pytest.mark.parametrize("name", ["lhaaso_1yr", "hawc_1yr", "veritas_50h", "hess_50h", "astri_50h"])
def test_instrument_sensitivities(name):
    """Sensitivity of other instruments: E² dN/dE between 1e-14 and 1e-10 TeV cm-2 s-1"""
    import astropy.units as u

    from sst1mpipe.performance.irf_plots import REFERENCE_DIRECTORY, read_irf_directory

    data = read_irf_directory(REFERENCE_DIRECTORY / name)
    assert set(data) == {"sensitivity"}
    sensitivity = data["sensitivity"]["ENERGY_FLUX_SENSITIVITY"].to_value(u.TeV / u.cm**2 / u.s)
    assert np.all((sensitivity > 1e-14) & (sensitivity < 1e-10))


@pytest.mark.parametrize("with_monitoring", [True, False])
def test_dl2_dl3(models, irfs, observation_monitoring, tmp_path, with_monitoring):
    """
    DL3 of an observation of telescope 21 (models and IRFs of telescope 1): read by gammapy with its
    point-like IRFs; IRFs of the closest zenith angle; pointing of the target without monitoring data
    """
    import shutil

    from astropy.coordinates import angular_separation
    from ctapipe.irf import OptimizationResult
    from gammapy.data import DataStore

    from sst1mpipe.io.dl3 import select_irf_directory
    from sst1mpipe.scripts.sst1mpipe_dl2_dl3_tool import DL2ToDL3

    run = OBSERVATION_RUNS[21]
    dl1 = tmp_path / "observation.dl1.h5"
    monitoring = monitoring_source_args(observation_monitoring[21]) if with_monitoring else []
    run_process(run["path"], dl1, f"--max-events={run['n_events']}", *monitoring, config=DATA_CONFIG_FILE)
    directory, _ = models
    dl2 = tmp_path / "observation.dl2.h5"
    reconstructors = [arg for path in sorted((directory / "nsb100").glob("*.pkl")) for arg in ("-r", str(path))]
    assert run_tool(ApplyModels(), argv=["-i", str(dl1), *reconstructors, "-o", str(dl2), "--no-progress"],
                    raises=True) == 0

    # IRFs at two zenith angles; the ze40 ones without gammaness cut (the test models are poor)
    irf, _ = irfs
    irf_directories = tmp_path / "irf"
    shutil.copytree(irf, irf_directories / "ze20")
    shutil.copytree(irf, irf_directories / "ze40")
    cuts = OptimizationResult.read(irf / "cuts.fits")
    cuts.gh_cuts["cut"] = 0.0
    cuts.write(irf_directories / "ze40" / "cuts.fits", overwrite=True)
    assert select_irf_directory(irf_directories, 30.5 * u.deg).name == "ze40"
    assert select_irf_directory(irf_directories, 25 * u.deg).name == "ze20"
    assert select_irf_directory(irf, 60 * u.deg) == irf

    output = tmp_path / "dl3"
    argv = ["-i", str(dl2), "--irf-directory", str(irf_directories / "ze40"), "-o", str(output)]
    assert run_tool(DL2ToDL3(), argv=argv, raises=True) == 0
    # GADF specification: the DL3 files and the index files
    from astropy.io import fits

    from sst1mpipe.io.dl3 import check_dl3_file, check_index

    dl3_files = sorted(output.glob("*.dl3.fits.gz"))
    assert len(dl3_files) == 1
    assert check_dl3_file(dl3_files[0]) == []
    assert check_index(output) == []
    header = fits.getheader(dl3_files[0], "EVENTS")
    assert header["INSTRUME"] == "SST-1M-1" and header["ORIGIN"] and header["CREATOR"].startswith("sst1mpipe-dl2-dl3")
    assert header["TIMESYS"] == "TT" and header["DATE-OBS"] and header["RA_OBJ"] == pytest.approx(header["RA_PNT"], abs=2)
    store = DataStore.from_dir(output)
    # the HDU types of GADF (the POINTING HDU is in the files, not in the index)
    assert set(store.hdu_table["HDU_TYPE"]) == {"events", "gti", "aeff", "edisp", "rad_max", "bkg"}
    assert {"OBS_ID", "OBS_MODE", "RA_PNT", "DEC_PNT", "TSTART", "TSTOP", "ZEN_PNT", "LIVETIME"} <= set(store.obs_table.colnames)
    observation = store.get_observations(required_irf="point-like")[0]
    events = observation.events.table
    valid = read_table(dl2, "/dl2/event/subarray/particle_type/RandomForestClassifier")["RandomForestClassifier_is_valid"]
    assert 0 < len(events) <= np.count_nonzero(valid)
    assert observation.obs_id == read_table(dl2, "/dl1/event/subarray/trigger")["obs_id"][0]
    assert observation.rad_max is not None and observation.aeff is not None
    # events around the pointing (wobble of the observation)
    pointing = observation.pointing.fixed_icrs
    separation = angular_separation(events["RA"].quantity, events["DEC"].quantity, pointing.ra, pointing.dec)
    assert np.all(separation < 5 * u.deg)
    target = read_table(dl2, "/configuration/observation/observation_block")["target"][0]
    assert store.obs_table["OBJECT"][0] == target
    assert 0 < store.obs_table["ZEN_PNT"].quantity[0] < 90 * u.deg

def test_cut_without_n_events():
    """the cuts of the optimization of the sensitivity have no n_events: all the bins drawn"""
    from astropy.table import QTable

    from sst1mpipe.performance.irf_plots import _cut

    table = QTable({"low": [1, 2] * u.TeV, "high": [2, 4] * u.TeV, "cut": [0.2, 0.1] * u.deg})
    np.testing.assert_allclose(_cut(table)[2], [0.2, 0.1])
    table["n_events"] = [5, 20]
    assert np.isnan(_cut(table)[2][0]) and _cut(table)[2][1] == 0.1


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

    # stereo direction of the Crab paper (StereoDispSignClosestCombiner of the stereo disp model): both telescopes
    disp = read_table(output, "/dl2/event/subarray/geometry/disp")
    n_telescopes = np.sum(np.asarray(disp["disp_telescopes"]), axis=1)
    assert np.any(disp["disp_is_valid"])
    assert np.all(n_telescopes[np.asarray(disp["disp_is_valid"], dtype=bool)] == 2)
    # the same models with another combiner (apply time): also the events with one valid telescope
    output_mean = tmp_path / "gamma_point_mean.dl2.h5"
    argv = ["-i", str(dl1_files["gamma_point"]), *reconstructors, "-o", str(output_mean), "--no-progress",
            "--ApplyModels.disp_stereo_combiner=StereoMeanCombiner"]
    assert run_tool(ApplyModels(), argv=argv, raises=True) == 0
    mean = read_table(output_mean, "/dl2/event/subarray/geometry/disp")
    assert np.count_nonzero(mean["disp_is_valid"]) >= np.count_nonzero(disp["disp_is_valid"])
