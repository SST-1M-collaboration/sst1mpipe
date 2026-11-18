"""
Plots of the models trained by the sst1mpipe training tools (`sst1mpipe.ml.tools`), made with
ctaplot: feature importances of the models, and performance from the cross validation of the
ctapipe training tools (``/cv_predictions/<telescope type>`` tables of the ``--cv-output`` file).
"""
import logging

import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from ctaplot import plots as ctaplot  # noqa: E402

__all__ = [
    "plot_feature_importances",
    "plot_energy_performance",
    "plot_classifier_performance",
    "plot_disp_performance",
]

logger = logging.getLogger(__name__)

#: true energy bins of the performance plots (TeV)
ENERGY_BINS = np.logspace(-1, 3, 21) * u.TeV


def _save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Plot saved in %s", path)
    return path


def _model_label(subarray, key):
    """Label of the model of a telescope type: its name and its telescopes"""
    tel_ids = subarray.get_tel_ids_for_type(key) if subarray is not None else []
    return f"{key} (tel {', '.join(map(str, tel_ids))})" if len(tel_ids) else str(key)


def plot_feature_importances(models, features, subarray, path, title):
    """
    Feature importances of the models of each telescope type (``models``: {key: scikit-learn model
    or tuple of models, e.g. the norm regressor and the sign classifier of the disp})

    Returns the path of the plot.
    """
    rows = []
    for key, model in models.items():
        parts = model if isinstance(model, tuple) else (model,)
        names = ("disp norm", "disp sign") if isinstance(model, tuple) else ("",)
        for name, part in zip(names, parts, strict=True):
            rows.append((f"{_model_label(subarray, key)} {name}".strip(), part.feature_importances_))
    fig, axes = plt.subplots(len(rows), 1, figsize=(max(8, 0.6 * len(features)), 3.5 * len(rows)), squeeze=False)
    for ax, (label, importances) in zip(axes[:, 0], rows, strict=True):
        ctaplot.plot_feature_importance(features, importances, ax=ax)
        ax.set_title(label)
        ax.set_ylabel("importance")
    fig.suptitle(title)
    fig.tight_layout()
    return _save(fig, path)


def _energy(table, column="true_energy"):
    """True energy of a cross validation table (TeV), None if not in it"""
    if column not in table.colnames:
        return None
    energy = table[column]
    return energy.quantity.to(u.TeV) if energy.unit is not None else energy * u.TeV


def plot_energy_performance(cv_tables, prefix, path):
    """Energy resolution, bias and migration matrix of the cross validation of the energy regressor"""
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for key, table in cv_tables.items():
        true = table["truth"].quantity.to(u.TeV) if table["truth"].unit else table["truth"] * u.TeV
        reco = table[f"{prefix}_energy"]
        reco = reco.quantity.to(u.TeV) if reco.unit else reco * u.TeV
        ctaplot.plot_energy_resolution(true, reco, ax=axes[0], bins=ENERGY_BINS, label=key)
        # ctaplot 0.6: energy_bias always uses its own bins (custom bins give a size mismatch)
        ctaplot.plot_energy_bias(true, reco, ax=axes[1], label=key)
        ctaplot.plot_migration_matrix(
            np.log10(true.to_value(u.TeV)), np.log10(reco.to_value(u.TeV)), ax=axes[2], colorbar=True,
            xy_line=True, hist2d_args=dict(bins=60, cmin=1),
        )
    axes[2].set(xlabel="log10(true energy / TeV)", ylabel="log10(reconstructed energy / TeV)", title="Migration matrix")
    for ax in axes[:2]:
        ax.legend()
    fig.suptitle(f"Energy regressor {prefix}: cross validation")
    fig.tight_layout()
    return _save(fig, path)


def plot_classifier_performance(cv_tables, prefix, path):
    """Gammaness distributions, ROC curves (all energies and per true energy) of the cross validation"""
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for key, table in cv_tables.items():
        truth = np.asarray(table["truth"])  # 1 for the positive class (gammas), 0 otherwise
        gammaness = np.asarray(table[f"{prefix}_prediction"])
        ctaplot.plot_gammaness_distribution(truth, gammaness, ax=axes[0], histtype="step", bins=50, density=True)
        ctaplot.plot_roc_curve_gammaness(truth, gammaness, gamma_label=1, ax=axes[1], label=key)
        energy = _energy(table)
        if energy is not None:
            ctaplot.plot_roc_curve_gammaness_per_energy(
                truth, gammaness, energy, gamma_label=1, energy_bins=ENERGY_BINS[::4], ax=axes[2],
            )
    axes[0].set(xlabel="gammaness", title="Gammaness (1: gammas, 0: background)")
    axes[1].legend()
    fig.suptitle(f"Particle classifier {prefix}: cross validation")
    fig.tight_layout()
    return _save(fig, path)


def plot_disp_performance(cv_tables, prefix, path):
    """Disp norm (migration, resolution per energy) and sign accuracy per energy of the cross validation"""
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for key, table in cv_tables.items():
        truth = table["truth"].quantity.to_value(u.deg) if table["truth"].unit else np.asarray(table["truth"])
        disp = table[f"{prefix}_parameter"]
        disp = disp.quantity.to_value(u.deg) if disp.unit else np.asarray(disp)
        ctaplot.plot_migration_matrix(
            np.abs(truth), np.abs(disp), ax=axes[0], colorbar=True, xy_line=True, hist2d_args=dict(bins=60, cmin=1),
        )
        energy = _energy(table)
        if energy is None:
            continue
        ctaplot.plot_resolution_per_energy(np.abs(truth), np.abs(disp), energy, ax=axes[1], bins=ENERGY_BINS, label=key)
        # fraction of the events with the right sign, per true energy
        index = np.digitize(energy.to_value(u.TeV), ENERGY_BINS.to_value(u.TeV)) - 1
        right = np.sign(truth) == np.sign(disp)
        centers = np.sqrt(ENERGY_BINS[:-1] * ENERGY_BINS[1:]).to_value(u.TeV)
        accuracy = [right[index == i].mean() if np.any(index == i) else np.nan for i in range(len(centers))]
        axes[2].plot(centers, accuracy, "o-", label=key)
    axes[0].set(xlabel="|true disp| (deg)", ylabel="|reconstructed disp| (deg)", title="Disp norm")
    axes[1].set(title="Resolution of the disp norm", ylabel="relative resolution of |disp| (68 %)")
    axes[2].set(xscale="log", xlabel="true energy (TeV)", ylabel="right sign fraction", ylim=(0, 1.05),
                title="Disp sign accuracy")
    for ax in axes[1:]:
        ax.legend()
    fig.suptitle(f"Disp reconstructor {prefix}: cross validation")
    fig.tight_layout()
    return _save(fig, path)
