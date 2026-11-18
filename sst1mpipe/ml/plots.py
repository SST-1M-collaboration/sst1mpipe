"""
Plots of the trained models (ctaplot): feature importances, and performance from the out-of-fold
predictions of the cross validation of the ctapipe training tools, per telescope id (one row each),
and in stereo (last row) for the array events with at least two telescopes, combined by the stereo
combiner of the model.
"""
import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.table import vstack  # noqa: E402
from ctapipe.reco.telescope_event_handling import get_subarray_index  # noqa: E402
from ctaplot import plots as ctaplot  # noqa: E402

from .stereo_combination import _candidates  # noqa: E402

__all__ = [
    "plot_feature_importances",
    "plot_energy_performance",
    "plot_classifier_performance",
    "plot_disp_performance",
]

#: true energy bins of the performance plots
ENERGY_BINS = np.logspace(-1, 3, 21) * u.TeV


def _quantity(column, unit):
    """A table column as a quantity (in ``unit`` if it has no unit)"""
    return column.quantity.to(unit) if column.unit is not None else u.Quantity(np.asarray(column), unit)


def _save(fig, path, title):
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _stereo(cv, combiner, columns, split=None):
    """
    Array events with at least two telescopes of the cross validation predictions ``cv``, combined
    by ``combiner``: ``columns`` {prediction column: telescope column of the combiner}, the stereo
    table has the prediction columns and the truth. ``split``: column of separate samples (the
    signal and background files of the classifier can have the same event ids).
    """
    samples = [cv] if split is None else [cv[cv[split] == value] for value in np.unique(cv[split])]
    stereo = []
    for table in samples:
        table = table.copy()
        table.sort(["obs_id", "event_id", "tel_id"])
        _, _, multiplicity, index = get_subarray_index(table)
        table = table[multiplicity[index] >= 2]
        if len(table) == 0:
            continue
        mono = table.copy()
        for name, tel_name in columns.items():
            mono[tel_name] = table[name]
        mono[f"{combiner.prefix}_tel_is_valid"] = np.ones(len(mono), dtype=bool)
        events = combiner.predict_table(mono)
        _, first = np.unique(get_subarray_index(mono)[3], return_index=True)
        for name in ("truth", "true_energy", "true_alt", "true_az"):
            if name in table.colnames:
                events[name] = table[name][first]
        stereo.append(events)
    return vstack(stereo) if stereo else None


def _selections(cv, stereo):
    """(label, table): one per telescope id, and the stereo events"""
    rows = [(f"tel {tel_id}", cv[cv["tel_id"] == tel_id]) for tel_id in np.unique(cv["tel_id"])]
    return rows + ([("stereo", stereo)] if stereo is not None and len(stereo) else [])


def _figure(n_rows, n_cols):
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.6 * n_cols, 4.3 * n_rows), squeeze=False)
    return fig, axes


def plot_feature_importances(models, features, path, title):
    """
    Feature importances of random forests, ``models``: {label: scikit-learn model, or (norm
    regressor, sign classifier) of the disp}
    """
    rows = []
    for label, model in models.items():
        if isinstance(model, tuple):
            rows += [(f"{label} disp norm", model[0]), (f"{label} disp sign", model[1])]
        else:
            rows.append((label, model))
    fig, axes = plt.subplots(len(rows), 1, figsize=(max(8, 0.6 * len(features)), 3.5 * len(rows)), squeeze=False)
    for ax, (label, model) in zip(axes[:, 0], rows, strict=True):
        ctaplot.plot_feature_importance(features, model.feature_importances_, ax=ax)
        ax.set(title=label, ylabel="importance")
    return _save(fig, path, title)


def plot_energy_performance(cv, model, path):
    """Energy resolution, bias and migration matrix"""
    column = f"{model.prefix}_energy"
    cv[column] = _quantity(cv[column], u.TeV)
    stereo = _stereo(cv, model.stereo_combiner, {column: f"{model.prefix}_tel_energy"})
    rows = _selections(cv, stereo)
    fig, axes = _figure(len(rows), 3)
    for (label, table), ax in zip(rows, axes, strict=True):
        true, reco = _quantity(table["truth"], u.TeV), _quantity(table[column], u.TeV)
        ctaplot.plot_energy_resolution(true, reco, ax=ax[0], bins=ENERGY_BINS)
        # ctaplot 0.6: energy_bias always uses its own bins (others give a size mismatch)
        ctaplot.plot_energy_bias(true, reco, ax=ax[1])
        ctaplot.plot_migration_matrix(
            np.log10(true.to_value(u.TeV)), np.log10(reco.to_value(u.TeV)), ax=ax[2], colorbar=True,
            xy_line=True, hist2d_args=dict(bins=60, cmin=1),
        )
        ax[0].set_title(f"{label}: energy resolution")
        ax[1].set_title(f"{label}: energy bias")
        ax[2].set(xlabel="log10(true energy / TeV)", ylabel="log10(reconstructed energy / TeV)",
                  title=f"{label}: migration matrix")
    return _save(fig, path, f"Energy regressor {model.prefix}: cross validation")


def plot_classifier_performance(cv, model, path):
    """Gammaness distributions and ROC curves (all energies, per true energy)"""
    column = f"{model.prefix}_prediction"
    stereo = _stereo(cv, model.stereo_combiner, {column: f"{model.prefix}_tel_prediction"}, split="truth")
    rows = _selections(cv, stereo)
    fig, axes = _figure(len(rows), 3)
    for (label, table), ax in zip(rows, axes, strict=True):
        truth, gammaness = np.asarray(table["truth"]), np.asarray(table[column])  # truth 1: gammas
        ctaplot.plot_gammaness_distribution(truth, gammaness, ax=ax[0], histtype="step", bins=50, density=True)
        ctaplot.plot_roc_curve_gammaness(truth, gammaness, gamma_label=1, ax=ax[1])
        ctaplot.plot_roc_curve_gammaness_per_energy(
            truth, gammaness, _quantity(table["true_energy"], u.TeV), gamma_label=1, energy_bins=ENERGY_BINS[::4],
            ax=ax[2],
        )
        ax[0].set(xlabel="gammaness", title=f"{label}: gammaness (1: gammas, 0: background)")
        ax[1].set_title(f"{label}: ROC curve")
        ax[2].set_title(f"{label}: ROC curve per true energy")
    return _save(fig, path, f"Particle classifier {model.prefix}: cross validation")


def plot_disp_performance(cv, model, path):
    """Disp norm migration and resolution, sign accuracy, angular resolution per true energy"""
    prefix = f"{model.prefix}"
    if "telescope_pointing_altitude" in cv.colnames:
        pointing = cv["telescope_pointing_altitude"].quantity, cv["telescope_pointing_azimuth"].quantity
    else:
        pointing = cv["subarray_pointing_lat"].quantity, cv["subarray_pointing_lon"].quantity
    cv[f"{prefix}_parameter"] = _quantity(cv[f"{prefix}_parameter"], u.deg)
    vectors = _candidates(
        cv[f"{prefix}_parameter"].quantity, cv["hillas_psi"].quantity, cv["hillas_fov_lon"].quantity,
        cv["hillas_fov_lat"].quantity, *pointing,
    )[:, 0]
    cv[f"{prefix}_alt"] = np.rad2deg(np.arcsin(vectors[:, 2])) * u.deg
    cv[f"{prefix}_az"] = np.rad2deg(np.arctan2(vectors[:, 1], vectors[:, 0]) % (2 * np.pi)) * u.deg
    stereo = _stereo(cv, model.stereo_combiner, {
        f"{prefix}_alt": f"{prefix}_tel_alt", f"{prefix}_az": f"{prefix}_tel_az",
        f"{prefix}_parameter": f"{prefix}_tel_parameter",
    })
    rows = _selections(cv, stereo)
    fig, axes = _figure(len(rows), 4)
    centers = np.sqrt(ENERGY_BINS[:-1] * ENERGY_BINS[1:]).to_value(u.TeV)
    for (label, table), ax in zip(rows, axes, strict=True):
        energy = _quantity(table["true_energy"], u.TeV)
        ctaplot.plot_angular_resolution_per_energy(
            _quantity(table["true_alt"], u.rad), _quantity(table[f"{prefix}_alt"], u.rad),
            _quantity(table["true_az"], u.rad), _quantity(table[f"{prefix}_az"], u.rad), energy,
            ax=ax[3], bins=ENERGY_BINS,
        )
        ax[3].set_title(f"{label}: angular resolution")
        if label == "stereo":
            for a in ax[:3]:
                a.axis("off")
            continue
        truth = _quantity(table["truth"], u.deg).to_value(u.deg)
        disp = table[f"{prefix}_parameter"].quantity.to_value(u.deg)
        ctaplot.plot_migration_matrix(
            np.abs(truth), np.abs(disp), ax=ax[0], colorbar=True, xy_line=True, hist2d_args=dict(bins=60, cmin=1),
        )
        ctaplot.plot_resolution_per_energy(np.abs(truth), np.abs(disp), energy, ax=ax[1], bins=ENERGY_BINS)
        index = np.digitize(energy.to_value(u.TeV), ENERGY_BINS.to_value(u.TeV)) - 1
        right = np.sign(truth) == np.sign(disp)
        ax[2].plot(centers, [right[index == i].mean() if np.any(index == i) else np.nan for i in range(len(centers))],
                   "o-")
        ax[0].set(xlabel="|true disp| (deg)", ylabel="|reconstructed disp| (deg)", title=f"{label}: disp norm")
        ax[1].set(title=f"{label}: disp norm resolution", ylabel="relative resolution of |disp| (68 %)")
        ax[2].set(xscale="log", xlabel="true energy (TeV)", ylabel="right sign fraction", ylim=(0, 1.05),
                  title=f"{label}: disp sign accuracy")
    return _save(fig, path, f"Disp reconstructor {prefix}: cross validation")
