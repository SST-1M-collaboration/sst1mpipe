"""
Training of the ctapipe reconstructors with the ctapipe tools (same options and configuration), and
plots of the trained models in ``plot_dir`` (`sst1mpipe.ml.plots`): feature importances, and
performance from the out-of-fold predictions of the cross validation, per telescope id and in stereo.

- ``sst1mpipe-train-energy-regressor``: ``ctapipe-train-energy-regressor``
- ``sst1mpipe-train-particle-classifier``: ``ctapipe-train-particle-classifier``
- ``sst1mpipe-train-disp-reconstructor``: ``ctapipe-train-disp-reconstructor``
"""
import sys
from pathlib import Path as PathlibPath
from unittest import mock

from astropy.table import hstack, vstack
from ctapipe.core import Tool
from ctapipe.core.traits import Path
from ctapipe.tools.train_disp_reconstructor import TrainDispReconstructor as CtapipeTrainDispReconstructor
from ctapipe.tools.train_energy_regressor import TrainEnergyRegressor as CtapipeTrainEnergyRegressor
from ctapipe.tools.train_particle_classifier import TrainParticleClassifier as CtapipeTrainParticleClassifier
from ctapipe.tools.utils import read_training_events

from sst1mpipe.ml import plots  # also makes the StereoMARSCombiner available

__all__ = ["TrainEnergyRegressor", "TrainParticleClassifier", "TrainDispReconstructor"]

ID_COLUMNS = ["obs_id", "event_id", "tel_id"]
#: columns of the weights of the stereo combiners
WEIGHT_COLUMNS = {
    "intensity": ["hillas_intensity"],
    "aspect-weighted-intensity": ["hillas_intensity", "hillas_length", "hillas_width"],
    "containment-weighted-intensity": [
        "hillas_intensity", "hillas_length", "hillas_width", "leakage_intensity_width_2",
    ],
}
#: columns of the training events kept with the cross validation predictions
PLOT_COLUMNS = ID_COLUMNS + sorted(set(sum(WEIGHT_COLUMNS.values(), []))) + [
    "true_energy", "true_alt", "true_az", "hillas_psi", "hillas_fov_lon", "hillas_fov_lat",
    "telescope_pointing_altitude", "telescope_pointing_azimuth", "subarray_pointing_lat", "subarray_pointing_lon",
]


class _PlotsMixin(Tool):
    """Plots of the models after the training of a ctapipe training tool (placed before it in the bases)"""

    plot_dir = Path(
        default_value=None,
        allow_none=True,
        file_ok=False,
        help="Directory of the plots of the models, by default <output name>_plots next to the output file",
    ).tag(config=True)

    #: attribute of the ctapipe tool with the trained reconstructor, and its performance plot
    model_attribute = None
    plot_performance = None

    def setup(self):
        output = PathlibPath(self.output_path)
        self.plot_dir = PathlibPath(self.plot_dir or output.parent / f"{output.stem}_plots")
        self.plot_dir.mkdir(parents=True, exist_ok=True)
        super().setup()
        self.model = getattr(self, self.model_attribute)
        self._training_events, self._cv_predictions = {}, []
        self._record_cross_validation()

    def start(self):
        # the ctapipe tools keep only the features: the ids are read too (and the stereo weights)
        ctapipe_tool = next(c for c in type(self).__mro__ if c.__module__.startswith("ctapipe.tools."))
        with mock.patch.object(sys.modules[ctapipe_tool.__module__], "read_training_events", self._read_events):
            super().start()

    def _read_events(self, *, feature_names, telescope_type, **kwargs):
        extra = ID_COLUMNS + WEIGHT_COLUMNS.get(self.model.stereo_combiner.weights, [])
        feature_names = feature_names + [c for c in extra if c not in feature_names]
        table = read_training_events(feature_names=feature_names, telescope_type=telescope_type, **kwargs)
        self._training_events[telescope_type] = table
        return table

    def _record_cross_validation(self):
        """Keeps the out-of-fold predictions of the CrossValidator with the ids of the telescope events"""
        cv = self.cross_validate
        split_data, fold_predictions = cv.split_data, cv.cross_validate
        test_indices = []

        def split(**kwargs):
            kfold = split_data(**kwargs)
            kfold_split = kfold.split

            def recording_split(*args):
                for train, test in kfold_split(*args):
                    test_indices.append(test)
                    yield train, test

            kfold.split = recording_split
            return kfold

        def predictions(telescope_type, train, test):
            result, metrics = fold_predictions(telescope_type, train, test)
            # the disp tool removes the ids from the table of the cross validation
            events = test if "tel_id" in test.colnames else self._training_events[telescope_type][test_indices[-1]]
            self._cv_predictions.append(hstack([events[[c for c in PLOT_COLUMNS if c in events.colnames]], result]))
            return result, metrics

        cv.split_data, cv.cross_validate = split, predictions

    def finish(self):
        super().finish()
        subarray = (self.signal_loader if hasattr(self, "signal_loader") else self.loader).subarray
        stem = PathlibPath(self.output_path).stem
        self.plots = [plots.plot_feature_importances(
            self.model._models, self.model.features, subarray, self.plot_dir / f"{stem}_feature_importance.png",
            title=f"{type(self.model).__name__}: feature importances",
        )]
        if self._cv_predictions:
            self.plots.append(self.plot_performance(
                vstack(self._cv_predictions), self.model, self.plot_dir / f"{stem}_performance.png",
            ))
        for path in self.plots:
            self.log.info("Plot saved in %s", path)


class TrainEnergyRegressor(_PlotsMixin, CtapipeTrainEnergyRegressor):
    """``ctapipe-train-energy-regressor`` with plots: feature importances, energy resolution, bias, migration"""

    name = "sst1mpipe-train-energy-regressor"
    aliases = {**CtapipeTrainEnergyRegressor.aliases, "plot-dir": "TrainEnergyRegressor.plot_dir"}
    model_attribute = "regressor"
    plot_performance = staticmethod(plots.plot_energy_performance)


class TrainParticleClassifier(_PlotsMixin, CtapipeTrainParticleClassifier):
    """``ctapipe-train-particle-classifier`` with plots: feature importances, gammaness, ROC curves"""

    name = "sst1mpipe-train-particle-classifier"
    aliases = {**CtapipeTrainParticleClassifier.aliases, "plot-dir": "TrainParticleClassifier.plot_dir"}
    model_attribute = "classifier"
    plot_performance = staticmethod(plots.plot_classifier_performance)


class TrainDispReconstructor(_PlotsMixin, CtapipeTrainDispReconstructor):
    """``ctapipe-train-disp-reconstructor`` with plots: feature importances, disp norm, sign accuracy"""

    name = "sst1mpipe-train-disp-reconstructor"
    aliases = {**CtapipeTrainDispReconstructor.aliases, "plot-dir": "TrainDispReconstructor.plot_dir"}
    model_attribute = "models"
    plot_performance = staticmethod(plots.plot_disp_performance)


def main_energy_regressor():
    TrainEnergyRegressor().run()


def main_particle_classifier():
    TrainParticleClassifier().run()


def main_disp_reconstructor():
    TrainDispReconstructor().run()
