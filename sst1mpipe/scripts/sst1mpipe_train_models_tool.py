"""
Training of the ctapipe reconstructors with the ctapipe tools (same options and configuration). At
the end, the plots of the trained models are made in ``plot_dir`` (`sst1mpipe.ml.plots`): feature
importances, and performance per telescope from the cross validation of ctapipe.

- ``sst1mpipe-train-energy-regressor``: ``ctapipe-train-energy-regressor``
- ``sst1mpipe-train-particle-classifier``: ``ctapipe-train-particle-classifier``
- ``sst1mpipe-train-disp-reconstructor``: ``ctapipe-train-disp-reconstructor``
"""
from pathlib import Path as PathlibPath

import numpy as np
from ctapipe.core.traits import Path
from ctapipe.io import read_table
from ctapipe.tools.train_disp_reconstructor import TrainDispReconstructor as CtapipeTrainDispReconstructor
from ctapipe.tools.train_energy_regressor import TrainEnergyRegressor as CtapipeTrainEnergyRegressor
from ctapipe.tools.train_particle_classifier import TrainParticleClassifier as CtapipeTrainParticleClassifier

from sst1mpipe.ml import plots

__all__ = ["TrainEnergyRegressor", "TrainParticleClassifier", "TrainDispReconstructor"]

PLOT_DIR_HELP = "Directory of the plots of the models, by default <output name>_plots next to the output file"


def setup_plots(tool):
    """
    Directory of the plots (``tool.plot_dir``, by default <output name>_plots next to the output
    file), and file of the cross validation in it if none is given (``--cv-output``)
    """
    output = PathlibPath(tool.output_path)
    tool.plot_dir = PathlibPath(tool.plot_dir or output.parent / f"{output.stem}_plots")
    tool.plot_dir.mkdir(parents=True, exist_ok=True)
    if tool.config.get("CrossValidator", {}).get("output_path") is None:
        tool.config.CrossValidator.output_path = str(tool.plot_dir / f"{output.stem}_cv.h5")


def cross_validation_per_telescope(cross_validator, subarray):
    """
    Predictions of the cross validation of ctapipe, ``{"tel<id>": table}``.

    ctapipe writes them per telescope type name (``/cv_predictions/<name>``): the SST-1M telescopes
    are two types of the same name (``SST_1M_DigiCam``), the folds of the second type follow the ones
    of the first one in the same table (``cv_fold`` back to 0).
    """
    h5file = getattr(cross_validator, "h5file", None)
    if h5file is None or "/cv_predictions" not in h5file:
        return {}
    h5file.flush()
    tables = {}
    for node in h5file.root.cv_predictions:
        table = read_table(h5file, node._v_pathname)
        types = [t for t in subarray.telescope_types if str(t) == node._v_name]
        folds = np.asarray(table["cv_fold"], dtype=np.int64)  # uint8 in the file
        starts = np.concatenate([[0], np.flatnonzero(np.diff(folds) < 0) + 1, [len(table)]])
        for tel_type, start, stop in zip(types, starts[:-1], starts[1:], strict=True):
            tel_ids = "_".join(map(str, subarray.get_tel_ids_for_type(tel_type)))
            tables[f"tel{tel_ids}"] = table[start:stop]
    return tables


def plot_feature_importances(tool, model, subarray, title):
    """Feature importances of the model of each telescope (type)"""
    models = {
        f"tel{'_'.join(map(str, subarray.get_tel_ids_for_type(tel_type)))}": tel_model
        for tel_type, tel_model in model._models.items()
    }
    return plots.plot_feature_importances(models, model.features, tool.plot_dir, PathlibPath(tool.output_path).stem,
                                          title=title)


class TrainEnergyRegressor(CtapipeTrainEnergyRegressor):
    """``ctapipe-train-energy-regressor`` with plots: feature importances, energy resolution, bias, migration"""

    name = "sst1mpipe-train-energy-regressor"
    plot_dir = Path(default_value=None, allow_none=True, file_ok=False, help=PLOT_DIR_HELP).tag(config=True)
    aliases = {**CtapipeTrainEnergyRegressor.aliases, "plot-dir": "TrainEnergyRegressor.plot_dir"}

    def setup(self):
        setup_plots(self)
        super().setup()

    def finish(self):
        subarray = self.loader.subarray
        self.plots = plot_feature_importances(self, self.regressor, subarray, "Energy regressor: feature importances")
        cv = cross_validation_per_telescope(self.cross_validate, subarray)
        if cv:
            self.plots += plots.plot_energy_performance(
                cv, self.regressor.prefix, self.plot_dir, PathlibPath(self.output_path).stem)
        super().finish()


class TrainParticleClassifier(CtapipeTrainParticleClassifier):
    """``ctapipe-train-particle-classifier`` with plots: feature importances, gammaness, ROC curves"""

    name = "sst1mpipe-train-particle-classifier"
    plot_dir = Path(default_value=None, allow_none=True, file_ok=False, help=PLOT_DIR_HELP).tag(config=True)
    aliases = {**CtapipeTrainParticleClassifier.aliases, "plot-dir": "TrainParticleClassifier.plot_dir"}

    def setup(self):
        setup_plots(self)
        super().setup()

    def finish(self):
        subarray = self.signal_loader.subarray
        self.plots = plot_feature_importances(self, self.classifier, subarray,
                                              "Particle classifier: feature importances")
        cv = cross_validation_per_telescope(self.cross_validate, subarray)
        if cv:
            self.plots += plots.plot_classifier_performance(
                cv, self.classifier.prefix, self.plot_dir, PathlibPath(self.output_path).stem)
        super().finish()


class TrainDispReconstructor(CtapipeTrainDispReconstructor):
    """``ctapipe-train-disp-reconstructor`` with plots: feature importances, disp norm, sign accuracy"""

    name = "sst1mpipe-train-disp-reconstructor"
    plot_dir = Path(default_value=None, allow_none=True, file_ok=False, help=PLOT_DIR_HELP).tag(config=True)
    aliases = {**CtapipeTrainDispReconstructor.aliases, "plot-dir": "TrainDispReconstructor.plot_dir"}

    def setup(self):
        setup_plots(self)
        super().setup()

    def finish(self):
        subarray = self.loader.subarray
        self.plots = plot_feature_importances(self, self.models, subarray, "Disp reconstructor: feature importances")
        cv = cross_validation_per_telescope(self.cross_validate, subarray)
        if cv:
            self.plots += plots.plot_disp_performance(
                cv, self.models.prefix, self.plot_dir, PathlibPath(self.output_path).stem)
        super().finish()


def main_energy_regressor():
    TrainEnergyRegressor().run()


def main_particle_classifier():
    TrainParticleClassifier().run()


def main_disp_reconstructor():
    TrainDispReconstructor().run()
