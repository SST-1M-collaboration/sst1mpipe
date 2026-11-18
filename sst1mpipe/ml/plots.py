"""
Plots of the trained models (ctaplot): feature importances, and performance from the cross validation
of the ctapipe training tools, one row per telescope.
"""
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

#: true energy bins of the performance plots
ENERGY_BINS = np.logspace(-1, 3, 21) * u.TeV


def _quantity(column, unit):
    """A table column as a quantity (in ``unit`` if it has no unit)"""
    return column.quantity.to(unit) if column.unit is not None else u.Quantity(np.asarray(column), unit)


def _figure(n_rows, n_cols):
    return plt.subplots(n_rows, n_cols, figsize=(5.6 * n_cols, 4.3 * n_rows), squeeze=False)


def _save(fig, path, title):
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


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


def plot_energy_performance(cv_tables, prefix, path):
    """Energy resolution, bias and migration matrix, ``cv_tables``: {label: cross validation table}"""
    fig, axes = _figure(len(cv_tables), 3)
    for (label, table), ax in zip(cv_tables.items(), axes, strict=True):
        true, reco = _quantity(table["truth"], u.TeV), _quantity(table[f"{prefix}_energy"], u.TeV)
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
    return _save(fig, path, f"Energy regressor {prefix}: cross validation")


def plot_classifier_performance(cv_tables, prefix, path):
    """Gammaness distributions and ROC curves (all energies, per true energy)"""
    fig, axes = _figure(len(cv_tables), 3)
    for (label, table), ax in zip(cv_tables.items(), axes, strict=True):
        truth, gammaness = np.asarray(table["truth"]), np.asarray(table[f"{prefix}_prediction"])  # truth 1: gammas
        ctaplot.plot_gammaness_distribution(truth, gammaness, ax=ax[0], histtype="step", bins=50, density=True)
        ctaplot.plot_roc_curve_gammaness(truth, gammaness, gamma_label=1, ax=ax[1])
        ctaplot.plot_roc_curve_gammaness_per_energy(
            truth, gammaness, _quantity(table["true_energy"], u.TeV), gamma_label=1, energy_bins=ENERGY_BINS[::4],
            ax=ax[2],
        )
        ax[0].set(xlabel="gammaness", title=f"{label}: gammaness (1: gammas, 0: background)")
        ax[1].set_title(f"{label}: ROC curve")
        ax[2].set_title(f"{label}: ROC curve per true energy")
    return _save(fig, path, f"Particle classifier {prefix}: cross validation")


def plot_disp_performance(cv_tables, prefix, path):
    """Disp norm migration and resolution, sign accuracy per true energy"""
    fig, axes = _figure(len(cv_tables), 3)
    centers = np.sqrt(ENERGY_BINS[:-1] * ENERGY_BINS[1:]).to_value(u.TeV)
    for (label, table), ax in zip(cv_tables.items(), axes, strict=True):
        truth = _quantity(table["truth"], u.deg).to_value(u.deg)
        disp = _quantity(table[f"{prefix}_parameter"], u.deg).to_value(u.deg)
        energy = _quantity(table["true_energy"], u.TeV)
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
