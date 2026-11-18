"""
sst1mpipe tools built on the ctapipe machine learning tools (``ctapipe.tools``), with the same
options and configuration sections (``TrainEnergyRegressor``, ``EnergyRegressor``, ...):

- `TrainEnergyRegressor`, `TrainParticleClassifier`, `TrainDispReconstructor`: after the training
  (``finish``), plots of the feature importances of the models and of their performance from the
  cross validation, made with ctaplot (`sst1mpipe.ml.plots`) in ``plot_dir``;
- `ApplyModels`: the models are chosen in a production of models with one subdirectory per node
  ``zeXX_azXX_nsbXX`` (zenith, azimuth, NSB rate): the node closest to the pointing of the input file.

The DL2 written is the one of ctapipe (column names, tables). Importing this module makes the
`StereoMARSCombiner` available (``stereo_combiner_cls`` of the reconstructors).
"""
import re
from pathlib import Path as PathlibPath

import astropy.units as u
import numpy as np
import tables
from ctapipe.core import Tool
from ctapipe.core.traits import Float, Path
from ctapipe.io import TableLoader, read_table
from ctapipe.tools.apply_models import ApplyModels as CtapipeApplyModels
from ctapipe.tools.train_disp_reconstructor import TrainDispReconstructor as CtapipeTrainDispReconstructor
from ctapipe.tools.train_energy_regressor import TrainEnergyRegressor as CtapipeTrainEnergyRegressor
from ctapipe.tools.train_particle_classifier import TrainParticleClassifier as CtapipeTrainParticleClassifier

from sst1mpipe.ml import plots
from sst1mpipe.ml.stereo_combination import StereoMARSCombiner  # noqa: F401 (registered for the reconstructors)

__all__ = [
    "TrainEnergyRegressor",
    "TrainParticleClassifier",
    "TrainDispReconstructor",
    "ApplyModels",
    "find_model_nodes",
    "closest_model_node",
    "main_train_energy_regressor",
    "main_train_particle_classifier",
    "main_train_disp_reconstructor",
    "main_apply_models",
]


class _PlotsMixin(Tool):
    """
    Plots of the trained models after the training of a ctapipe training tool (a Tool, for its
    trait ``plot_dir``: placed before the ctapipe tool in the bases of the sst1mpipe tools)
    """

    plot_dir = Path(
        default_value=None,
        allow_none=True,
        file_ok=False,
        help=(
            "Directory of the plots of the models (feature importances, performance from the cross"
            " validation). By default <output name>_plots next to the output file."
        ),
    ).tag(config=True)

    def _setup_plots(self):
        if self.output_path is None:
            return
        output = PathlibPath(self.output_path)
        if self.plot_dir is None:
            self.plot_dir = output.parent / f"{output.stem}_plots"
        PathlibPath(self.plot_dir).mkdir(parents=True, exist_ok=True)
        # the performance plots need the predictions of the cross validation: written in the plot
        # directory if no --cv-output is given
        cross_validator = self.config.get("CrossValidator", {})
        if cross_validator.get("output_path") is None and cross_validator.get("n_cross_validations", 5) > 0:
            self.config.CrossValidator.output_path = str(PathlibPath(self.plot_dir) / f"{output.stem}_cv.h5")

    def setup(self):
        self._setup_plots()
        super().setup()

    def finish(self):
        super().finish()
        self.plots = self.make_plots()
        for path in self.plots:
            self.log.info("Plot saved in %s", path)

    def _cv_tables(self):
        """Predictions of the cross validation ``{telescope type: table}`` (empty without cross validation)"""
        path = self.cross_validate.output_path
        if path is None or not PathlibPath(path).exists():
            return {}
        with tables.open_file(path) as h5file:
            names = [node._v_name for node in h5file.walk_nodes("/cv_predictions", "Table")] \
                if "/cv_predictions" in h5file else []
        return {name: read_table(path, f"/cv_predictions/{name}") for name in names}

    def _plot_path(self, name):
        return PathlibPath(self.plot_dir) / f"{PathlibPath(self.output_path).stem}_{name}.png"

    def make_plots(self):
        """Feature importances and performance plots, returns their paths"""
        component = self.model_component
        paths = [plots.plot_feature_importances(
            component._models, component.features, self.loader.subarray, self._plot_path("feature_importance"),
            title=f"{type(component).__name__}: feature importances",
        )]
        cv_tables = self._cv_tables()
        if cv_tables:
            paths.append(self.plot_performance(cv_tables, component.prefix, self._plot_path("performance")))
        return paths


class TrainEnergyRegressor(_PlotsMixin, CtapipeTrainEnergyRegressor):
    """``ctapipe-train-energy-regressor``, with plots of the feature importances and of the energy resolution, bias and migration"""

    name = "sst1mpipe-train-energy-regressor"
    aliases = {**CtapipeTrainEnergyRegressor.aliases, "plot-dir": "TrainEnergyRegressor.plot_dir"}

    @property
    def model_component(self):
        return self.regressor

    @staticmethod
    def plot_performance(cv_tables, prefix, path):
        return plots.plot_energy_performance(cv_tables, prefix, path)


class TrainParticleClassifier(_PlotsMixin, CtapipeTrainParticleClassifier):
    """``ctapipe-train-particle-classifier``, with plots of the feature importances, gammaness and ROC curves"""

    name = "sst1mpipe-train-particle-classifier"
    aliases = {**CtapipeTrainParticleClassifier.aliases, "plot-dir": "TrainParticleClassifier.plot_dir"}

    @property
    def model_component(self):
        return self.classifier

    @property
    def loader(self):
        return self.signal_loader

    @staticmethod
    def plot_performance(cv_tables, prefix, path):
        return plots.plot_classifier_performance(cv_tables, prefix, path)


class TrainDispReconstructor(_PlotsMixin, CtapipeTrainDispReconstructor):
    """``ctapipe-train-disp-reconstructor``, with plots of the feature importances, disp norm and sign accuracy"""

    name = "sst1mpipe-train-disp-reconstructor"
    aliases = {**CtapipeTrainDispReconstructor.aliases, "plot-dir": "TrainDispReconstructor.plot_dir"}

    @property
    def model_component(self):
        return self.models

    @staticmethod
    def plot_performance(cv_tables, prefix, path):
        return plots.plot_disp_performance(cv_tables, prefix, path)


#: name of the subdirectories of a production of models: zenith, azimuth (deg), NSB rate (MHz)
MODEL_NODE_PATTERN = re.compile(r"^ze(?P<zenith>[-\d.]+)_az(?P<azimuth>[-\d.]+)_nsb(?P<nsb>[-\d.]+)$")


def find_model_nodes(models_directory):
    """
    The nodes of a production of models: ``{subdirectory: (zenith, azimuth, nsb)}`` of the
    subdirectories named ``zeXX_azXX_nsbXX`` (the other subdirectories are ignored)
    """
    nodes = {}
    for path in sorted(PathlibPath(models_directory).iterdir()):
        match = MODEL_NODE_PATTERN.match(path.name)
        if path.is_dir() and match:
            nodes[path] = tuple(float(match[key]) for key in ("zenith", "azimuth", "nsb"))
    return nodes


def closest_model_node(nodes, zenith, azimuth=None, nsb_rate=None):
    """
    The node of a production of models (`find_model_nodes`) closest to a pointing: the closest NSB
    rate if ``nsb_rate`` is given (the models are made for given NSB levels), then the closest
    zenith, then the closest azimuth if ``azimuth`` is given (angles in deg, NSB in MHz)
    """
    if len(nodes) == 0:
        raise ValueError("No model node (subdirectory zeXX_azXX_nsbXX)")
    candidates = list(nodes.items())
    if nsb_rate is not None:
        best = min(abs(nsb - nsb_rate) for _, (_, _, nsb) in candidates)
        candidates = [c for c in candidates if abs(c[1][2] - nsb_rate) == best]
    best = min(abs(ze - zenith) for _, (ze, _, _) in candidates)
    candidates = [c for c in candidates if abs(c[1][0] - zenith) == best]
    if azimuth is not None:
        def azimuth_distance(candidate):
            return abs((candidate[1][1] - azimuth + 180) % 360 - 180)
        candidates.sort(key=azimuth_distance)
    return candidates[0][0]


class ApplyModels(CtapipeApplyModels):
    """
    ``ctapipe-apply-models`` with the models of the node of a production of models closest to the
    pointing of the input file: ``models_directory`` has one subdirectory per node,
    ``zeXX_azXX_nsbXX`` (zenith and azimuth in deg, NSB rate in MHz), with the trained models
    (``*.pkl``) of the node. The node is the closest in NSB rate (if ``nsb_rate`` is given), then in
    zenith, then in azimuth, to the median pointing of the telescope events of the input file.
    With ``--reconstructor`` (``reconstructor_paths``), these models are used instead.

    The output is the DL2 of ctapipe. The `StereoMARSCombiner` can be used by the models.
    """

    name = "sst1mpipe-apply-models"

    models_directory = Path(
        default_value=None,
        allow_none=True,
        exists=True,
        file_ok=False,
        help="Production of models: one subdirectory zeXX_azXX_nsbXX per node, with its models (*.pkl)",
    ).tag(config=True)

    nsb_rate = Float(
        default_value=None,
        allow_none=True,
        help="NSB rate of the input file (MHz) to choose the models; the NSB is not used if not given",
    ).tag(config=True)

    aliases = {
        **CtapipeApplyModels.aliases,
        ("m", "models-directory"): "ApplyModels.models_directory",
        "nsb-rate": "ApplyModels.nsb_rate",
    }

    def setup(self):
        if len(self.reconstructor_paths) == 0 and self.models_directory is not None and self.input_url is not None:
            node = self.closest_node()
            self.reconstructor_paths = sorted(node.glob("*.pkl"))
            self.log.info("Models of the node %s: %s", node.name, [p.name for p in self.reconstructor_paths])
            if len(self.reconstructor_paths) == 0:
                self.log.critical("No model (*.pkl) in %s", node)
                self.exit(1)
        super().setup()

    def closest_node(self):
        """The node of ``models_directory`` closest to the pointing of the input file"""
        zenith, azimuth = input_pointing(self.input_url)
        self.log.info("Pointing of the input file: zenith %.2f deg, azimuth %.2f deg", zenith, azimuth)
        return closest_model_node(find_model_nodes(self.models_directory), zenith, azimuth, self.nsb_rate)


def input_pointing(path):
    """Median pointing (zenith, azimuth in deg) of the telescope events of a DL1 file"""
    with TableLoader(path) as loader:
        events = loader.read_telescope_events(dl1_parameters=False, pointing=True, observation_info=True)
    if "telescope_pointing_altitude" in events.colnames:
        altitude, azimuth = events["telescope_pointing_altitude"], events["telescope_pointing_azimuth"]
    else:
        altitude, azimuth = events["subarray_pointing_lat"], events["subarray_pointing_lon"]
    altitude = u.Quantity(altitude, u.deg).to_value(u.deg)
    azimuth = u.Quantity(azimuth, u.deg).to_value(u.deg)
    # median azimuth around the circular mean (pointings near azimuth 0 / 360)
    mean_azimuth = np.rad2deg(np.arctan2(np.nanmean(np.sin(np.deg2rad(azimuth))), np.nanmean(np.cos(np.deg2rad(azimuth)))))
    azimuth = mean_azimuth + np.nanmedian((azimuth - mean_azimuth + 180) % 360 - 180)
    return 90 - float(np.nanmedian(altitude)), float(azimuth % 360)


def main_train_energy_regressor():
    TrainEnergyRegressor().run()


def main_train_particle_classifier():
    TrainParticleClassifier().run()


def main_train_disp_reconstructor():
    TrainDispReconstructor().run()


def main_apply_models():
    ApplyModels().run()
