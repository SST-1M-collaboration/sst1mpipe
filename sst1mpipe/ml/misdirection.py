"""
Misdirection regressor and event classes (placeholder, work in progress).

Plan:

1. ``MisdirectionRegressor``: a scikit-learn regression model per telescope type, as the ctapipe
   reconstructors (``SKLearnRegressionReconstructor``), of the angular distance between the true and
   the reconstructed directions (target ``log10(misdirection / deg)``), from a set of features: DL1
   parameters, disp norm and sign score, reconstructed energy, ``disp_goodness_of_fit`` of the MARS
   combination... Trained on diffuse gammas (DL2 of ``sst1mpipe-apply-models``) by a training tool
   as ``sst1mpipe-train-*``, applied by ``sst1mpipe-apply-models``. The prediction could fill the
   uncertainty of the geometry of ctapipe (``disp_ang_distance_uncert``) instead of a new column.
2. Event classes: binning of the predicted misdirection (``misdirection_edges``: quantiles of the
   point gammas, per reconstructed energy bin), ``event_type`` of each event.
3. DL3: ``EVENT_TYPE`` in the EVENTS table, and IRFs per event class.
"""
__all__ = ["MisdirectionRegressor", "misdirection_edges", "event_classes"]


class MisdirectionRegressor:
    """
    TODO: ``SKLearnRegressionReconstructor`` with ``target = "misdirection"`` (log10 of the angular
    distance between the true and the reconstructed directions), ``predict_table`` and the event
    loop ``__call__``, and its training tool with plots (``sst1mpipe.ml.plots``).
    """

    def __init__(self, *args, **kwargs):
        raise NotImplementedError("misdirection regressor: work in progress")


def misdirection_edges(misdirection, reco_energy, energy_bins, quantiles=(0.25, 0.5, 0.75)):
    """
    TODO: edges of the event classes, per reconstructed energy bin: quantiles of the predicted
    misdirection of the point gammas, shape (n_energy_bins, len(quantiles))
    """
    raise NotImplementedError("event classes: work in progress")


def event_classes(misdirection, reco_energy, energy_bins, edges):
    """
    TODO: event class (0: best direction) of each event: bin of its predicted misdirection with the
    edges of its reconstructed energy bin, -1 outside the energy bins
    """
    raise NotImplementedError("event classes: work in progress")

