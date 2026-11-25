"""Training of the ctapipe reconstructors with the sst1mpipe tools (mono and stereo models)"""
import pytest
from ctapipe.core import run_tool
from ctapipe.reco import Reconstructor

from conftest import run_process
from sst1mpipe.resources import MC_CONFIG_FILES, STEREO_TRAIN_CONFIG_FILES, TEST_DATA_DIR, TRAIN_CONFIG_FILES
from sst1mpipe.scripts.sst1mpipe_train_models_tool import (
    TrainDispReconstructor,
    TrainEnergyRegressor,
    TrainParticleClassifier,
)

SIMTEL_DIR = TEST_DATA_DIR / "simtel"
GAMMA_DIFFUSE = SIMTEL_DIR / "gamma_200_800E3GeV_20_20deg_ATM52_102168.corsika.gz.NSBmed4.simtel.gz"
GAMMA_POINT = SIMTEL_DIR / "gamma_point_200_800E3GeV_20_20deg_ATM52_300669.corsika.gz.NSBmed4.simtel.gz"
PROTON = SIMTEL_DIR / "proton_400_1300E3GeV_20_20deg_ATM52_204508.corsika.gz.NSBmed4.simtel.gz"


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


def test_stereo_models(dl1_files, tmp_path):
    """
    stereo models: the stereo parameters of the HillasReconstructor (impact distance, h_max) as features,
    trained on the events reconstructed by both telescopes
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
