"""
Plots of the trained models: feature importances, and performance from the cross validation of the
ctapipe training tools. One plot per figure, per telescope (``<stem>_<label>_<quantity>.png``), with
the same axis limits for all the figures (``LIMITS``) to compare them.
"""
import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from ctaplot import plots as ctaplot  # noqa: E402

__all__ = [
    "LIMITS",
    "ENERGY_BINS",
    "bias_resolution",
    "figure",
    "save",
    "plot_feature_importances",
    "plot_energy_performance",
    "plot_classifier_performance",
    "plot_disp_performance",
]

#: true energy bins of the performance plots (5 per decade)
ENERGY_BINS = np.logspace(-1, np.log10(300), 18) * u.TeV

#: axis limits of all the figures
LIMITS = {
    "energy": (0.1, 300),  # TeV
    "bias": (-0.5, 0.5),  # relative
    "resolution": (0, 0.5),  # relative
    "migration": (0.1, 300),  # TeV, reconstructed vs true energy
    "disp": (0, 4),  # deg
    "fraction": (0, 1.02),
}

INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
plt.rcParams.update({
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.titlecolor": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2, "lines.markersize": 5,
    "font.size": 11, "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def _quantity(column, unit):
    """A table column as a quantity (in ``unit`` if it has no unit)"""
    return column.quantity.to(unit) if column.unit is not None else u.Quantity(np.asarray(column), unit)


def figure():
    """One plot per figure"""
    return plt.subplots(figsize=(7, 5))


def save(fig, path, title):
    fig.axes[0].set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def bias_resolution(true, reco, energy, bins=ENERGY_BINS, min_events=10):
    """
    Bias (median) and resolution (half of the 16-84 % interval) of the relative error
    ``reco / true - 1`` per energy bin, as the benchmarks of pyirf / ctapipe; NaN with fewer than
    ``min_events`` events
    """
    error = np.asarray(reco / true - 1, dtype=float)
    index = np.digitize(energy.to_value(u.TeV), bins.to_value(u.TeV)) - 1
    bias, resolution = np.full(len(bins) - 1, np.nan), np.full(len(bins) - 1, np.nan)
    for i in range(len(bins) - 1):
        values = error[(index == i) & np.isfinite(error)]
        if len(values) >= min_events:
            q16, q50, q84 = np.percentile(values, [16, 50, 84])
            bias[i], resolution[i] = q50, (q84 - q16) / 2
    return bias, resolution


def _per_energy(ax, values, bins=ENERGY_BINS, **kwargs):
    """Values per energy bin: points with the bin widths"""
    edges = bins.to_value(u.TeV)
    centers = np.sqrt(edges[:-1] * edges[1:])
    ok = np.isfinite(values)
    ax.errorbar(centers[ok], values[ok], xerr=[(centers - edges[:-1])[ok], (edges[1:] - centers)[ok]], fmt="o",
                elinewidth=2, capsize=0, **kwargs)
    ax.set(xscale="log", xlim=LIMITS["energy"])


def plot_feature_importances(models, features, directory, stem, title):
    """
    Feature importances of random forests, one figure per model: ``models`` {label: scikit-learn
    model, or (norm regressor, sign classifier) of the disp}
    """
    paths = []
    for label, model in models.items():
        parts = {"": model} if not isinstance(model, tuple) else {"_norm": model[0], "_sign": model[1]}
        for suffix, part in parts.items():
            fig, ax = plt.subplots(figsize=(max(7, 0.55 * len(features)), 5))
            ctaplot.plot_feature_importance(features, part.feature_importances_, ax=ax)
            ax.set(ylabel="importance", ylim=(0, 1))
            name = {"": "", "_norm": " (disp norm)", "_sign": " (disp sign)"}[suffix]
            paths.append(save(fig, directory / f"{stem}_{label}_feature_importance{suffix}.png",
                              f"{title}, {label}{name}"))
    return paths


def plot_energy_performance(cv_tables, prefix, directory, stem):
    """Energy resolution, bias and migration matrix, ``cv_tables``: {label: cross validation table}"""
    paths = []
    for label, table in cv_tables.items():
        true, reco = _quantity(table["truth"], u.TeV), _quantity(table[f"{prefix}_energy"], u.TeV)
        bias, resolution = bias_resolution(true, reco, true)

        fig, ax = figure()
        _per_energy(ax, resolution, color=COLORS[0])
        ax.set(xlabel="true energy (TeV)", ylabel="energy resolution", ylim=LIMITS["resolution"])
        paths.append(save(fig, directory / f"{stem}_{label}_resolution.png", f"Energy resolution, {label}"))

        fig, ax = figure()
        _per_energy(ax, bias, color=COLORS[0])
        ax.axhline(0, color=MUTED, lw=1)
        ax.set(xlabel="true energy (TeV)", ylabel="energy bias", ylim=LIMITS["bias"])
        paths.append(save(fig, directory / f"{stem}_{label}_bias.png", f"Energy bias, {label}"))

        fig, ax = figure()
        edges = np.geomspace(*LIMITS["migration"], 61)
        counts, _, _ = np.histogram2d(true.to_value(u.TeV), reco.to_value(u.TeV), bins=[edges, edges])
        mesh = ax.pcolormesh(edges, edges, np.ma.masked_equal(counts.T, 0), cmap="Blues")
        ax.plot(LIMITS["migration"], LIMITS["migration"], color=MUTED, lw=1, ls="--")
        ax.set(xscale="log", yscale="log", xlim=LIMITS["migration"], ylim=LIMITS["migration"],
               xlabel="true energy (TeV)", ylabel="reconstructed energy (TeV)")
        fig.colorbar(mesh, ax=ax, label="events")
        paths.append(save(fig, directory / f"{stem}_{label}_migration.png", f"Energy migration, {label}"))
    return paths


def plot_classifier_performance(cv_tables, prefix, directory, stem):
    """Gammaness distributions and ROC curves (all energies, per true energy)"""
    paths = []
    for label, table in cv_tables.items():
        truth, gammaness = np.asarray(table["truth"]), np.asarray(table[f"{prefix}_prediction"])  # truth 1: gammas

        fig, ax = figure()
        bins = np.linspace(0, 1, 41)
        for value, name, color in ((1, "gammas", COLORS[0]), (0, "background", COLORS[1])):
            ax.hist(gammaness[truth == value], bins=bins, histtype="step", density=True, lw=2, color=color, label=name)
        ax.set(xlabel="gammaness", ylabel="probability density", xlim=(0, 1))
        ax.legend(frameon=False)
        paths.append(save(fig, directory / f"{stem}_{label}_gammaness.png", f"Gammaness, {label}"))

        fig, ax = figure()
        ctaplot.plot_roc_curve_gammaness(truth, gammaness, gamma_label=1, ax=ax, color=COLORS[0])
        ax.set(xlim=LIMITS["fraction"], ylim=LIMITS["fraction"])
        paths.append(save(fig, directory / f"{stem}_{label}_roc.png", f"ROC curve, {label}"))

        fig, ax = figure()
        ctaplot.plot_roc_curve_gammaness_per_energy(
            truth, gammaness, _quantity(table["true_energy"], u.TeV), gamma_label=1,
            energy_bins=ENERGY_BINS[::4], ax=ax,
        )
        ax.set(xlim=LIMITS["fraction"], ylim=LIMITS["fraction"])
        paths.append(save(fig, directory / f"{stem}_{label}_roc_per_energy.png", f"ROC curves per true energy, {label}"))
    return paths


def plot_disp_performance(cv_tables, prefix, directory, stem):
    """Disp norm migration, bias and resolution of the disp norm, sign accuracy per true energy"""
    paths = []
    for label, table in cv_tables.items():
        truth = _quantity(table["truth"], u.deg).to_value(u.deg)
        disp = _quantity(table[f"{prefix}_parameter"], u.deg).to_value(u.deg)
        energy = _quantity(table["true_energy"], u.TeV)

        fig, ax = figure()
        edges = np.linspace(*LIMITS["disp"], 61)
        counts, _, _ = np.histogram2d(np.abs(truth), np.abs(disp), bins=[edges, edges])
        mesh = ax.pcolormesh(edges, edges, np.ma.masked_equal(counts.T, 0), cmap="Blues")
        ax.plot(LIMITS["disp"], LIMITS["disp"], color=MUTED, lw=1, ls="--")
        ax.set(xlim=LIMITS["disp"], ylim=LIMITS["disp"], xlabel="|true disp| (deg)", ylabel="|reconstructed disp| (deg)")
        fig.colorbar(mesh, ax=ax, label="events")
        paths.append(save(fig, directory / f"{stem}_{label}_disp_migration.png", f"Disp norm, {label}"))

        bias, resolution = bias_resolution(np.abs(truth), np.abs(disp), energy)
        fig, ax = figure()
        _per_energy(ax, resolution, color=COLORS[0])
        ax.set(xlabel="true energy (TeV)", ylabel="disp norm resolution", ylim=LIMITS["resolution"])
        paths.append(save(fig, directory / f"{stem}_{label}_disp_resolution.png", f"Disp norm resolution, {label}"))

        fig, ax = figure()
        _per_energy(ax, bias, color=COLORS[0])
        ax.axhline(0, color=MUTED, lw=1)
        ax.set(xlabel="true energy (TeV)", ylabel="disp norm bias", ylim=LIMITS["bias"])
        paths.append(save(fig, directory / f"{stem}_{label}_disp_bias.png", f"Disp norm bias, {label}"))

        index = np.digitize(energy.to_value(u.TeV), ENERGY_BINS.to_value(u.TeV)) - 1
        right = np.sign(truth) == np.sign(disp)
        accuracy = np.array([right[index == i].mean() if np.count_nonzero(index == i) >= 10 else np.nan
                             for i in range(len(ENERGY_BINS) - 1)])
        fig, ax = figure()
        _per_energy(ax, accuracy, color=COLORS[0])
        ax.set(xlabel="true energy (TeV)", ylabel="right sign fraction", ylim=LIMITS["fraction"])
        paths.append(save(fig, directory / f"{stem}_{label}_disp_sign.png", f"Disp sign accuracy, {label}"))
    return paths
