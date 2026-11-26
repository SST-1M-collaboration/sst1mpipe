"""
Plots of the IRFs and benchmarks of ``ctapipe-compute-irf`` (point-like), and of the gamma/hadron
separation of the DL2: one plot per figure, with the axis limits of `IRF_LIMITS` for all the figures.
Each function takes ``inputs`` {label: data} and draws one series per label: one label for the plots
of one mode (e.g. ``tel1``), several to compare them. ``references``: other IRF directories of the same
format drawn as lines, e.g. the published performance in `REFERENCE_DIRECTORY` (only the quantities
published, see `write_reference_directory`):

- ``crab_paper_2025/<mode>/ze<zenith>``: SST-1M Crab paper (arXiv:2506.01733), tel1 and stereo
  at ze30, ze50 (energy bias and resolution, angular resolution, sensitivity), tel1, tel2 and
  stereo at ze20 (effective area)
- ``lhaaso_1yr``, ``hawc_1yr``, ``veritas_50h``, ``hess_50h``, ``astri_50h``: sensitivity of other
  instruments
"""
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import QTable

from pyirf.spectral import CRAB_MAGIC_JHEAP2015

from sst1mpipe.ml.plots import COLORS, LIMITS, MUTED, figure, save

__all__ = ["IRF_LIMITS", "REFERENCE_DIRECTORY", "read_irf_directory", "write_reference_directory", "plot_irfs",
           "plot_gamma_hadron_separation"]

#: published performance in IRF directories (scripts/reference)
REFERENCE_DIRECTORY = Path(__file__).parents[1] / "resources" / "reference"

#: axis limits of the IRF figures, in addition to the ones of the models (`sst1mpipe.ml.plots.LIMITS`)
IRF_LIMITS = {
    **LIMITS,
    "energy": (0.1, 1e3),  # TeV, true and reconstructed energy
    "bias": (-20, 20),  # %
    "resolution": (0, 30),  # %
    "migration": (0, 2),  # reconstructed / true energy
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
    The point-like IRFs, benchmarks and cuts of ``ctapipe-compute-irf`` and
    ``ctapipe-optimize-event-selection`` in a directory (``irf_point_like.fits.gz``, ``benchmarks_point_like.fits.gz``, ``cuts.fits``): the quantities
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


#: axis labels of the IRF figures (matplotlib mathtext: the exponents do not depend on the glyphs of the font)
E_TRUE = r"$E_\mathrm{true}$ (TeV)"
E_RECO = r"$E_\mathrm{reco}$ (TeV)"
LABELS = {
    "effective_area": r"$A_\mathrm{eff}$ (m$^2$)",
    "effective_area_protons": r"$A_\mathrm{eff}$, protons (m$^2$)",
    "energy_bias": "energy bias (%)",
    "energy_resolution": "energy resolution (%)",
    "angular_resolution": r"$\theta_{68}$ (deg)",
    "theta_cut": r"$\theta$ cut (deg)",
    "gammaness_cut": "gammaness cut",
    "background": r"background (s$^{-1}$ sr$^{-1}$ TeV$^{-1}$)",
    "sensitivity_crab": "sensitivity (% Crab)",
    "sensitivity": r"$E^2\,\mathrm{d}N/\mathrm{d}E$ (TeV cm$^{-2}$ s$^{-1}$)",
    "migration": r"$E_\mathrm{reco}/E_\mathrm{true}$",
}


def _reference(ax, low, high, values, color, label):
    """Values per energy bin of a reference: steps (gaps without value)"""
    x = np.ravel(np.column_stack([low, high]))
    y = np.repeat(np.where(np.isfinite(values) & (values != 0), values, np.nan), 2)
    ax.plot(x, y, color=color, lw=1.5, ls="--", label=label)


def _column(name, values, unit=None):
    """Column of one row of the HDUs of ctapipe-compute-irf: values (n_energy) or (offset, n_energy)"""
    values = np.atleast_2d(values) if np.ndim(values) > 1 else np.asarray(values, dtype=float)
    dim = "(" + ",".join(str(n) for n in values.shape[::-1]) + ")"
    return fits.Column(name, format=f"{values.size}D", dim=dim, unit=unit, array=values[np.newaxis])


def _hdu(name, bins, columns, reference, comment):
    columns = [_column("ENERG_LO", bins[:-1], "TeV"), _column("ENERG_HI", bins[1:], "TeV"),
               _column("THETA_LO", [0.0], "deg"), _column("THETA_HI", [1.0], "deg"), *columns]
    hdu = fits.BinTableHDU.from_columns(columns, name=name)
    hdu.header["REFERENC"] = reference
    if comment:
        hdu.header["COMMENT"] = comment
    return hdu


def write_reference_directory(directory, reference, true_energy_bins=None, reco_energy_bins=None,
                              effective_area=None, energy_bias=None, energy_resolution=None,
                              angular_resolution=None, sensitivity=None, aeff_energy_bins=None, comments=None):
    """
    Published performance in an IRF directory of the format of ``ctapipe-compute-irf`` (point-like),
    with the quantities given only (NaN: bins without value), for `plot_irfs` (``references``)

    Parameters
    ----------
    directory: path
    reference: str, publication (header keyword REFERENC)
    true_energy_bins, reco_energy_bins, aeff_energy_bins: energy bins (TeV) of the benchmarks per true
        energy (bias, resolutions), per reconstructed energy (sensitivity), of the effective area
    effective_area: m2 per ``aeff_energy_bins`` bin
    energy_bias, energy_resolution, angular_resolution (deg): per ``true_energy_bins`` bin
    sensitivity: E² dN/dE (TeV cm-2 s-1) per ``reco_energy_bins`` bin (relative: to the Crab of MAGIC, `_crab`)
    comments: {quantity: comment} (definitions of the publication)
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    comments = comments or {}

    def row(values):
        return np.asarray(values, dtype=float)[np.newaxis]

    if effective_area is not None:
        hdus = [_hdu("EFFECTIVE AREA", np.asarray(aeff_energy_bins), [_column("EFFAREA", row(effective_area), "m2")],
                     reference, comments.get("effective_area", ""))]
        fits.HDUList([fits.PrimaryHDU(), *hdus]).writeto(directory / "irf_point_like.fits.gz", overwrite=True)
    hdus = []
    if energy_bias is not None or energy_resolution is not None:
        n = len(true_energy_bins) - 1
        bias = np.full(n, np.nan) if energy_bias is None else energy_bias
        resolution = np.full(n, np.nan) if energy_resolution is None else energy_resolution
        hdus.append(_hdu("ENERGY BIAS RESOLUTION", np.asarray(true_energy_bins),
                         [_column("BIAS", row(bias)), _column("RESOLUTION", row(resolution))], reference,
                         comments.get("energy", "true energy")))
    if angular_resolution is not None:
        hdus.append(_hdu("ANGULAR RESOLUTION", np.asarray(true_energy_bins),
                         [_column("ANGULAR_RESOLUTION_68", row(angular_resolution), "deg")], reference,
                         comments.get("angular_resolution", "true energy")))
    if sensitivity is not None:
        bins = np.asarray(reco_energy_bins)
        energy = np.sqrt(bins[:-1] * bins[1:])
        sensitivity = np.asarray(sensitivity, dtype=float)
        energy_flux = (sensitivity * u.Unit("TeV cm-2 s-1")).to_value("erg s-1 cm-2")
        flux = (sensitivity / energy**2 * u.Unit("TeV-1 cm-2 s-1")).to_value("erg-1 s-1 cm-2")
        hdus.append(_hdu("SENSITIVITY", bins, [
            _column("RELATIVE_SENSITIVITY", row(sensitivity / _crab(energy))),
            _column("FLUX_SENSITIVITY", row(flux), "erg-1 s-1 cm-2"),
            _column("ENERGY_FLUX_SENSITIVITY", row(energy_flux), "erg s-1 cm-2"),
        ], reference, comments.get("sensitivity", "reconstructed energy")))
    if hdus:
        fits.HDUList([fits.PrimaryHDU(), *hdus]).writeto(directory / "benchmarks_point_like.fits.gz", overwrite=True)


def _cut(table, min_events=10):
    """
    Cut per reconstructed energy bin, NaN for the bins with the fill value (too few events: the percentile
    cuts have the column ``n_events``, not the cuts of the optimization of the sensitivity)
    """
    values = np.asarray(u.Quantity(table["cut"]).value, dtype=float)
    if "n_events" in table.colnames:
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
    """E² dN/dE of the Crab Nebula (MAGIC, Aleksić et al. 2015, as the SST-1M Crab paper and the sst1mpipe IRF tools) in TeV cm-2 s-1"""
    energy = u.Quantity(energy, u.TeV)
    return (energy**2 * CRAB_MAGIC_JHEAP2015(energy)).to_value(u.TeV / u.cm**2 / u.s)


def plot_irfs(inputs, directory, title="", references=None):
    """
    One figure per quantity in ``directory`` (``<quantity>.png``), one series per input:
    ``inputs`` {label: data of `read_irf_directory`}, ``references`` (same format) drawn as lines. The
    quantities of none of the inputs are not drawn. No figure titles: the labels of the series and
    ``title`` (e.g. ``30°, optimized cuts``, as the title of the legend) say what is drawn.
    """
    references = references or {}
    colors = dict(zip([*inputs, *(f"reference {label}" for label in references)], COLORS * 2, strict=False))
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = []

    def legend(ax):
        if ax.get_legend_handles_labels()[0]:
            ax.legend(frameon=False, fontsize=8, title=title or None, title_fontsize=8)

    def draw(ax, key, get):
        """Series of the inputs and references with ``key``: False if none"""
        drawn = False
        for label, data in inputs.items():
            if key in data:
                _series(ax, *get(data), colors[label], label)
                drawn = True
        for label, data in references.items():
            if key in data:
                _reference(ax, *get(data), colors[f"reference {label}"], label)
                drawn = True
        return drawn

    def per_energy(name, xlabel, ylim, get, key, log=False, zero_line=False):
        fig, ax = figure()
        if not draw(ax, key, get):
            plt.close(fig)
            return
        if zero_line:
            ax.axhline(0, color=MUTED, lw=1)
        ax.set(xscale="log", xlim=IRF_LIMITS["energy"], ylim=ylim, xlabel=xlabel, ylabel=LABELS[name])
        if log:
            ax.set_yscale("log")
        legend(ax)
        paths.append(save(fig, directory / f"{name}.png"))

    per_energy("effective_area", E_TRUE, IRF_LIMITS["effective_area"],
               lambda d: (*_edges(d["aeff"]), d["aeff"]["EFFAREA"].ravel().to_value(u.m**2)), "aeff", log=True)
    per_energy("effective_area_protons", E_TRUE, IRF_LIMITS["effective_area"],
               lambda d: (*_edges(d["aeff_protons"]), d["aeff_protons"]["EFFAREA"].ravel().to_value(u.m**2)),
               "aeff_protons", log=True)
    # bias and resolution in %
    per_energy("energy_bias", E_TRUE, IRF_LIMITS["bias"],
               lambda d: (*_edges(d["energy"]), 100 * np.ravel(d["energy"]["BIAS"])), "energy", zero_line=True)
    per_energy("energy_resolution", E_TRUE, IRF_LIMITS["resolution"],
               lambda d: (*_edges(d["energy"]), 100 * np.ravel(d["energy"]["RESOLUTION"])), "energy")
    per_energy("angular_resolution", E_TRUE, IRF_LIMITS["angular_resolution"],
               lambda d: (*_edges(d["angular"]), d["angular"]["ANGULAR_RESOLUTION_68"].ravel().to_value(u.deg)),
               "angular")
    # cuts: bins with at least 10 events (the other ones have the fill value of ctapipe)
    per_energy("theta_cut", E_RECO, IRF_LIMITS["theta_cut"], lambda d: _cut(d["theta_cuts"]), "theta_cuts")
    per_energy("gammaness_cut", E_RECO, IRF_LIMITS["gammaness_cut"], lambda d: _cut(d["gh_cuts"]), "gh_cuts")
    per_energy("background", E_RECO, IRF_LIMITS["background"],
               lambda d: (*_edges(d["background"]),
                          d["background"]["BKG"].ravel().to_value(u.Unit("s-1 sr-1 TeV-1"))),
               "background", log=True)
    per_energy("sensitivity_crab", E_RECO, IRF_LIMITS["crab"],
               lambda d: (*_edges(d["sensitivity"]), 100 * np.ravel(d["sensitivity"]["RELATIVE_SENSITIVITY"])),
               "sensitivity", log=True)

    # E² flux sensitivity, with the Crab Nebula
    fig, ax = figure()
    energy = np.geomspace(*IRF_LIMITS["energy"], 50)
    for fraction in (1, 0.1, 0.01):
        ax.plot(energy, fraction * _crab(energy), color=MUTED, lw=1, ls=":")
        ax.annotate(f"{fraction:.0%} Crab", (0.15, fraction * _crab(0.15)), color=MUTED, fontsize=9,
                    xytext=(0, 4), textcoords="offset points")
    if draw(ax, "sensitivity", lambda d: (*_edges(d["sensitivity"]), d["sensitivity"]["ENERGY_FLUX_SENSITIVITY"]
                                          .ravel().to_value(u.TeV / u.cm**2 / u.s))):
        ax.set(xscale="log", yscale="log", xlim=IRF_LIMITS["energy"], ylim=IRF_LIMITS["sensitivity"], xlabel=E_RECO,
               ylabel=LABELS["sensitivity"])
        legend(ax)
        paths.append(save(fig, directory / "sensitivity.png"))
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
        ax.set(xscale="log", xlim=IRF_LIMITS["energy"], ylim=IRF_LIMITS["migration"], xlabel=E_TRUE,
               ylabel=LABELS["migration"])
        fig.colorbar(mesh, ax=ax, label="probability density")
        ax.text(0.03, 0.95, ", ".join(text for text in (label, title) if text), transform=ax.transAxes, fontsize=8,
                va="top")
        name = "energy_migration" if len(inputs) == 1 else f"energy_migration_{label}"
        paths.append(save(fig, directory / f"{name}.png"))
    return paths


def plot_gamma_hadron_separation(inputs, directory, prefix="RandomForestClassifier", title=""):
    """
    Gammaness of the gammas and protons, and ROC curve (gamma: the positive class; true positive rate:
    gammas classified as gammas, false positive rate: protons classified as gammas) with its area
    (AUC), one series per input: ``inputs``
    {label: (gamma DL2 file, proton DL2 file)}
    """
    from ctapipe.io import read_table
    from sklearn.metrics import roc_auc_score, roc_curve

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
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
        ax.legend(frameon=False, fontsize=8, title=", ".join(text for text in (label, title) if text), title_fontsize=8)
        name = "gammaness" if len(inputs) == 1 else f"gammaness_{label}"
        paths.append(save(fig, directory / f"{name}.png"))

    fig, ax = figure()
    for (label, (gammas, protons)), color in zip(gammaness.items(), COLORS, strict=False):
        # gamma: the positive class
        truth = np.concatenate([np.ones(len(gammas)), np.zeros(len(protons))])
        scores = np.concatenate([gammas, protons])
        false_positive, true_positive, _ = roc_curve(truth, scores)
        ax.plot(false_positive, true_positive, color=color, label=f"{label} (AUC {roc_auc_score(truth, scores):.3f})")
    ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls="--")
    ax.set(xlim=IRF_LIMITS["fraction"], ylim=IRF_LIMITS["fraction"], xlabel="false positive rate (protons)",
           ylabel="true positive rate (gammas)")
    ax.legend(frameon=False, fontsize=8, title=title or None, title_fontsize=8)
    paths.append(save(fig, directory / "roc.png"))
    return paths
