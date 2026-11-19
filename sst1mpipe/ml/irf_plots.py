"""
Plots of the IRFs and benchmarks of ``ctapipe-compute-irf`` (point-like), and of the gamma/hadron
separation of the DL2: one plot per figure, with the axis limits of `IRF_LIMITS` for all the figures.
Each function takes ``inputs`` {label: data} and draws one series per label: one label for the plots
of one mode (e.g. ``tel1``), several to compare them. ``references``: other IRF directories of the same
format drawn as lines, e.g. the performance of the SST-1M Crab paper (arXiv:2506.01733) in
`REFERENCE_DIRECTORY` (``<mode>/ze<zenith>``: only the quantities of the paper).
"""
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import QTable

from .plots import COLORS, LIMITS, MUTED, figure, save

__all__ = ["IRF_LIMITS", "REFERENCE_DIRECTORY", "read_irf_directory", "plot_irfs", "plot_gamma_hadron_separation"]

#: performance of the SST-1M Crab paper (Alispach et al. 2025, arXiv:2506.01733) extracted from its
#: figures (scripts/reference/extract_crab_paper_irfs.py), in IRF directories <mode>/ze<zenith>:
#: tel1, strict_stereo (ze30, ze50: benchmarks), tel1, tel2, strict_stereo (ze20: effective area)
REFERENCE_DIRECTORY = Path(__file__).parents[1] / "resources" / "reference" / "crab_paper_2025"

#: axis limits of the IRF figures, in addition to the ones of the models (`sst1mpipe.ml.plots.LIMITS`)
IRF_LIMITS = {
    **LIMITS,
    "effective_area": (1e0, 1e7),  # m2
    "angular_resolution": (0, 0.5),  # deg
    "theta_cut": (0, 0.6),  # deg
    "gammaness_cut": (0, 1),
    "background": (1e-6, 1e2),  # s-1 sr-1 TeV-1
    "sensitivity": (1e-13, 1e-10),  # TeV cm-2 s-1
    "crab": (1, 1e3),  # % Crab
}


def read_irf_directory(directory):
    """
    The point-like IRFs, benchmarks and cuts of a directory of ``mc_irf.bash``
    (``irf_point_like.fits.gz``, ``benchmarks_point_like.fits.gz``, ``cuts.fits``): the quantities
    of the files and HDUs found (e.g. a directory of `REFERENCE_DIRECTORY`: only some of them)
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"No IRF directory {directory}")
    data = {}
    files = {
        "irf_point_like.fits.gz": (("EFFECTIVE AREA", "aeff"), ("EFFECTIVE AREA PROTONS", "aeff_protons"),
                                   ("ENERGY DISPERSION", "edisp"), ("RAD_MAX", "rad_max"),
                                   ("BACKGROUND", "background")),
        "benchmarks_point_like.fits.gz": (("ENERGY BIAS RESOLUTION", "energy"), ("ANGULAR RESOLUTION", "angular"),
                                          ("SENSITIVITY", "sensitivity")),
        "cuts.fits": (("GH_CUTS", "gh_cuts"), ("RAD_MAX", "theta_cuts")),
    }
    for filename, hdus in files.items():
        if not (directory / filename).is_file():
            continue
        with fits.open(directory / filename) as hdul:
            for hdu, name in hdus:
                if hdu in hdul:
                    table = QTable.read(hdul[hdu], mask_invalid=False)
                    data[name] = table if filename == "cuts.fits" else table[0]
    return data


def _reference(ax, low, high, values, color, label):
    """Values per energy bin of a reference: steps (gaps without value)"""
    x = np.ravel(np.column_stack([low, high]))
    y = np.repeat(np.where(np.isfinite(values) & (values != 0), values, np.nan), 2)
    ax.plot(x, y, color=color, lw=1.5, ls="--", label=label)


def _cut(table, min_events=10):
    """Cut per reconstructed energy bin, NaN for the bins with the fill value (too few events)"""
    values = np.asarray(u.Quantity(table["cut"]).value, dtype=float)
    values[np.asarray(table["n_events"]) < min_events] = np.nan
    return table["low"].to_value(u.TeV), table["high"].to_value(u.TeV), values


def _edges(row, prefix="ENERG"):
    return row[f"{prefix}_LO"].ravel().to_value(u.TeV), row[f"{prefix}_HI"].ravel().to_value(u.TeV)


def _series(ax, low, high, values, color, label):
    """Values per energy bin: points with the bin widths"""
    centers = np.sqrt(low * high)
    ok = np.isfinite(values) & (values != 0)
    ax.errorbar(centers[ok], values[ok], xerr=[(centers - low)[ok], (high - centers)[ok]], fmt="o", color=color,
                elinewidth=2, capsize=0, label=label)


def _crab(energy):
    """E² dN/dE of the Crab Nebula (HEGRA, Aharonian et al. 2004) in TeV cm-2 s-1"""
    return (energy**2 * 2.83e-11 * energy**-2.62)


def plot_irfs(inputs, directory, title="", references=None):
    """
    One figure per quantity in ``directory`` (``<quantity>.png``), one series per input:
    ``inputs`` {label: data of `read_irf_directory`}, ``references`` (same format) drawn as lines. The
    quantities of none of the inputs are not drawn.
    """
    references = references or {}
    colors = dict(zip([*inputs, *(f"reference {label}" for label in references)], COLORS * 2, strict=False))
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    suffix = f", {title}" if title else ""
    paths = []
    labels = len(inputs) + len(references) > 1

    def legend(ax):
        if ax.get_legend_handles_labels()[0]:
            ax.legend(frameon=False, fontsize=8)

    def draw(ax, key, get):
        """Series of the inputs and references with ``key``: False if none"""
        drawn = False
        for label, data in inputs.items():
            if key in data:
                _series(ax, *get(data), colors[label], label if labels else None)
                drawn = True
        for label, data in references.items():
            if key in data:
                _reference(ax, *get(data), colors[f"reference {label}"], label)
                drawn = True
        return drawn

    def per_energy(name, ylabel, xlabel, ylim, get, key, log=False, zero_line=False, plot_title=None):
        fig, ax = figure()
        if not draw(ax, key, get):
            plt.close(fig)
            return
        if zero_line:
            ax.axhline(0, color=MUTED, lw=1)
        ax.set(xscale="log", xlim=IRF_LIMITS["energy"], ylim=ylim, xlabel=xlabel, ylabel=ylabel)
        if log:
            ax.set_yscale("log")
        legend(ax)
        paths.append(save(fig, directory / f"{name}.png", (plot_title or ylabel.capitalize()) + suffix))

    true, reco = "true energy (TeV)", "reconstructed energy (TeV)"
    per_energy("effective_area", "effective area (m²)", true, IRF_LIMITS["effective_area"],
               lambda d: (*_edges(d["aeff"]), d["aeff"]["EFFAREA"].ravel().to_value(u.m**2)), "aeff", log=True,
               plot_title="Effective area (gammas)")
    per_energy("effective_area_protons", "effective area (m²)", true, IRF_LIMITS["effective_area"],
               lambda d: (*_edges(d["aeff_protons"]), d["aeff_protons"]["EFFAREA"].ravel().to_value(u.m**2)),
               "aeff_protons", log=True, plot_title="Effective area (protons)")
    per_energy("energy_bias", "energy bias", true, IRF_LIMITS["bias"],
               lambda d: (*_edges(d["energy"]), np.ravel(d["energy"]["BIAS"])), "energy", zero_line=True)
    per_energy("energy_resolution", "energy resolution", true, IRF_LIMITS["resolution"],
               lambda d: (*_edges(d["energy"]), np.ravel(d["energy"]["RESOLUTION"])), "energy")
    per_energy("angular_resolution", "68 % containment radius (deg)", true, IRF_LIMITS["angular_resolution"],
               lambda d: (*_edges(d["angular"]), d["angular"]["ANGULAR_RESOLUTION_68"].ravel().to_value(u.deg)),
               "angular", plot_title="Angular resolution")
    # cuts: bins with at least 10 events (the other ones have the fill value of ctapipe)
    per_energy("theta_cut", "θ cut (deg)", reco, IRF_LIMITS["theta_cut"], lambda d: _cut(d["theta_cuts"]),
               "theta_cuts", plot_title="Direction cut")
    per_energy("gammaness_cut", "gammaness cut", reco, IRF_LIMITS["gammaness_cut"], lambda d: _cut(d["gh_cuts"]), "gh_cuts")
    per_energy("background", "rate (s⁻¹ sr⁻¹ TeV⁻¹)", reco, IRF_LIMITS["background"],
               lambda d: (*_edges(d["background"]),
                          d["background"]["BKG"].ravel().to_value(u.Unit("s-1 sr-1 TeV-1"))),
               "background", log=True, plot_title="Background rate")
    per_energy("sensitivity_crab", "flux sensitivity (% Crab)", reco, IRF_LIMITS["crab"],
               lambda d: (*_edges(d["sensitivity"]), 100 * np.ravel(d["sensitivity"]["RELATIVE_SENSITIVITY"])),
               "sensitivity", log=True, plot_title="Sensitivity, 50 h")

    # E² flux sensitivity, with the Crab Nebula
    fig, ax = figure()
    energy = np.geomspace(*IRF_LIMITS["energy"], 50)
    for fraction in (1, 0.1, 0.01):
        ax.plot(energy, fraction * _crab(energy), color=MUTED, lw=1, ls=":")
        ax.annotate(f"{fraction:.0%} Crab", (0.15, fraction * _crab(0.15)), color=MUTED, fontsize=9,
                    xytext=(0, 4), textcoords="offset points")
    if draw(ax, "sensitivity", lambda d: (*_edges(d["sensitivity"]), d["sensitivity"]["ENERGY_FLUX_SENSITIVITY"]
                                          .ravel().to_value(u.TeV / u.cm**2 / u.s))):
        ax.set(xscale="log", yscale="log", xlim=IRF_LIMITS["energy"], ylim=IRF_LIMITS["sensitivity"], xlabel=reco,
               ylabel="E² dN/dE (TeV cm⁻² s⁻¹)")
        legend(ax)
        paths.append(save(fig, directory / "sensitivity.png", f"Flux sensitivity, 50 h{suffix}"))
    else:
        plt.close(fig)

    # energy migration: one figure per input
    for label, data in inputs.items():
        if "edisp" not in data:
            continue
        row = data["edisp"]
        low, high = _edges(row)
        migra = np.concatenate([row["MIGRA_LO"].ravel(), row["MIGRA_HI"].ravel()[-1:]])
        matrix = np.asarray(row["MATRIX"])
        matrix = matrix.reshape(-1, len(migra) - 1, len(low))[0]  # (offset, migra, energy): first offset bin
        fig, ax = figure()
        mesh = ax.pcolormesh(np.concatenate([low, high[-1:]]), migra, np.ma.masked_equal(matrix, 0), cmap="Blues")
        ax.axhline(1, color=MUTED, lw=1, ls="--")
        ax.set(xscale="log", xlim=IRF_LIMITS["energy"], ylim=(migra[0], migra[-1]), xlabel=true,
               ylabel="reconstructed / true energy")
        fig.colorbar(mesh, ax=ax, label="probability density")
        name = "energy_migration" if len(inputs) == 1 else f"energy_migration_{label}"
        paths.append(save(fig, directory / f"{name}.png", f"Energy migration, {label}{suffix}"))
    return paths


def plot_gamma_hadron_separation(inputs, directory, prefix="RandomForestClassifier", title=""):
    """
    Gammaness of the gammas and protons, and ROC curve, one series per input: ``inputs``
    {label: (gamma DL2 file, proton DL2 file)}
    """
    from ctapipe.io import read_table
    from sklearn.metrics import roc_curve

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    suffix = f", {title}" if title else ""
    gammaness = {}
    for label, (gamma_file, proton_file) in inputs.items():
        values = []
        for path in (gamma_file, proton_file):
            table = read_table(path, f"/dl2/event/subarray/particle_type/{prefix}")
            valid = np.asarray(table[f"{prefix}_is_valid"], dtype=bool)
            values.append(np.asarray(table[f"{prefix}_prediction"])[valid])
        gammaness[label] = values

    paths = []
    bins = np.linspace(0, 1, 41)
    for label, (gammas, protons) in gammaness.items():
        fig, ax = figure()
        for values, name, color in ((gammas, "gammas", COLORS[0]), (protons, "protons", COLORS[1])):
            ax.hist(values, bins=bins, histtype="step", density=True, lw=2, color=color, label=name)
        ax.set(xlim=(0, 1), xlabel="gammaness", ylabel="probability density")
        ax.legend(frameon=False)
        name = "gammaness" if len(inputs) == 1 else f"gammaness_{label}"
        paths.append(save(fig, directory / f"{name}.png", f"Gammaness, {label}{suffix}"))

    fig, ax = figure()
    for (label, (gammas, protons)), color in zip(gammaness.items(), COLORS, strict=False):
        truth = np.concatenate([np.ones(len(gammas)), np.zeros(len(protons))])
        false_positive, true_positive, _ = roc_curve(truth, np.concatenate([gammas, protons]))
        ax.plot(false_positive, true_positive, color=color, label=label)
    ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls="--")
    ax.set(xlim=IRF_LIMITS["fraction"], ylim=IRF_LIMITS["fraction"], xlabel="proton efficiency",
           ylabel="gamma efficiency")
    ax.legend(frameon=False, fontsize=8)
    paths.append(save(fig, directory / "roc.png", f"ROC curve{suffix}"))
    return paths
