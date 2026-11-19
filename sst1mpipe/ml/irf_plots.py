"""
Plots of the IRFs and benchmarks of ``ctapipe-compute-irf`` (point-like), and of the gamma/hadron
separation of the DL2: one plot per figure, with the axis limits of `IRF_LIMITS` for all the figures.
Each function takes ``inputs`` {label: data} and draws one series per label: one label for the plots
of one mode (e.g. ``tel1``), several to compare them. ``reference``: published performance drawn as
lines with the series (`read_reference`, default: the SST-1M Crab paper, arXiv:2506.01733).
"""
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import QTable, Table

from .plots import COLORS, LIMITS, MUTED, figure, save

__all__ = ["IRF_LIMITS", "REFERENCE_FILE", "read_irf_directory", "read_reference", "plot_irfs",
           "plot_gamma_hadron_separation"]

#: performance of the SST-1M Crab paper (Alispach et al. 2025, arXiv:2506.01733), extracted from its
#: figures (scripts/reference/extract_crab_paper_irfs.py)
REFERENCE_FILE = Path(__file__).parents[1] / "resources" / "reference" / "sst1m_crab_paper_2025.ecsv"
#: reference lines: color per regime, line style per zenith angle
REFERENCE_COLORS = {"SST-1M-1 mono": "#6b6b6b", "SST-1M-2 mono": "#a8a8a8", "stereo": "#1a1a1a"}
REFERENCE_STYLES = {20: "-", 30: "--", 50: ":"}

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


def _row(path, hdu):
    return QTable.read(path, hdu=hdu)[0]


def read_irf_directory(directory):
    """
    The point-like IRFs, benchmarks and cuts of a directory of ``mc_irf.bash``
    (``irf_point_like.fits.gz``, ``benchmarks_point_like.fits.gz``, ``cuts.fits``)
    """
    directory = Path(directory)
    irf, benchmarks = directory / "irf_point_like.fits.gz", directory / "benchmarks_point_like.fits.gz"
    data = {}
    for hdu, name in (("EFFECTIVE AREA", "aeff"), ("EFFECTIVE AREA PROTONS", "aeff_protons"),
                      ("ENERGY DISPERSION", "edisp"), ("RAD_MAX", "rad_max"), ("BACKGROUND", "background")):
        data[name] = _row(irf, hdu)
    for hdu, name in (("ENERGY BIAS RESOLUTION", "energy"), ("ANGULAR RESOLUTION", "angular"),
                      ("SENSITIVITY", "sensitivity")):
        data[name] = _row(benchmarks, hdu)
    with fits.open(directory / "cuts.fits") as cuts:
        data["gh_cuts"] = QTable.read(cuts["GH_CUTS"])
        data["theta_cuts"] = QTable.read(cuts["RAD_MAX"])
    return data


def read_reference(path=None, zeniths=None, stage="analysis level"):
    """
    Reference performance (default: `REFERENCE_FILE`): table with the columns ``quantity``,
    ``series``, ``regime``, ``zenith``, ``stage``, ``x_low``, ``x_high`` (equal for the curves),
    ``value``. ``zeniths``: only these zenith angles (deg) if given; ``stage``: the effective area at
    this stage of the analysis only (paper: triggered, survived cleaning, quality cuts, analysis level)
    """
    table = Table.read(path or REFERENCE_FILE, format="ascii.ecsv")
    keep = np.ones(len(table), dtype=bool)
    if zeniths:
        keep &= np.isin(table["zenith"], list(zeniths))
    if stage:
        keep &= (table["quantity"] != "effective_area") | (table["stage"] == stage)
    return table[keep]


def _reference(ax, reference, quantity, transform=None):
    """Lines of the reference series of ``quantity``: steps for the binned ones"""
    if reference is None:
        return
    table = reference[reference["quantity"] == quantity]
    for series in dict.fromkeys(table["series"]):
        rows = table[table["series"] == series]
        low, high, values = (np.asarray(rows[name], dtype=float) for name in ("x_low", "x_high", "value"))
        if transform is not None:
            values = transform(low, high, values)
        style = {"color": REFERENCE_COLORS.get(rows["regime"][0], MUTED),
                 "ls": REFERENCE_STYLES.get(int(rows["zenith"][0]), "-."), "lw": 1.5, "label": f"paper: {series}"}
        if np.all(low == high):
            ax.plot(low, values, **style)
        else:
            # contiguous bins (NaN between separate groups)
            x, y = [], []
            for i in range(len(rows)):
                if i and not np.isclose(low[i], high[i - 1]):
                    x.append(np.nan), y.append(np.nan)
                x += [low[i], high[i]]
                y += [values[i], values[i]]
            ax.plot(x, y, **style)


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


def plot_irfs(inputs, directory, title="", reference=None):
    """
    One figure per quantity in ``directory`` (``<quantity>.png``), one series per input:
    ``inputs`` {label: data of `read_irf_directory`}, ``reference``: table of `read_reference`
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    suffix = f", {title}" if title else ""
    paths = []
    labels = len(inputs) > 1 or reference is not None

    def legend(ax):
        if ax.get_legend_handles_labels()[0]:
            ax.legend(frameon=False, fontsize=8)

    def per_energy(name, ylabel, xlabel, ylim, get, log=False, zero_line=False, plot_title=None, quantity=None,
                   transform=None):
        fig, ax = figure()
        for (label, data), color in zip(inputs.items(), COLORS, strict=False):
            low, high, values = get(data)
            _series(ax, low, high, values, color, label if labels else None)
        _reference(ax, reference, quantity or name, transform)
        if zero_line:
            ax.axhline(0, color=MUTED, lw=1)
        ax.set(xscale="log", xlim=IRF_LIMITS["energy"], ylim=ylim, xlabel=xlabel, ylabel=ylabel)
        if log:
            ax.set_yscale("log")
        legend(ax)
        paths.append(save(fig, directory / f"{name}.png", (plot_title or ylabel.capitalize()) + suffix))

    true, reco = "true energy (TeV)", "reconstructed energy (TeV)"
    per_energy("effective_area", "effective area (m²)", true, IRF_LIMITS["effective_area"],
               lambda d: (*_edges(d["aeff"]), d["aeff"]["EFFAREA"].ravel().to_value(u.m**2)), log=True,
               plot_title="Effective area (gammas)")
    per_energy("effective_area_protons", "effective area (m²)", true, IRF_LIMITS["effective_area"],
               lambda d: (*_edges(d["aeff_protons"]), d["aeff_protons"]["EFFAREA"].ravel().to_value(u.m**2)),
               log=True, plot_title="Effective area (protons)")
    per_energy("energy_bias", "energy bias", true, IRF_LIMITS["bias"],
               lambda d: (*_edges(d["energy"]), np.ravel(d["energy"]["BIAS"])), zero_line=True)
    per_energy("energy_resolution", "energy resolution", true, IRF_LIMITS["resolution"],
               lambda d: (*_edges(d["energy"]), np.ravel(d["energy"]["RESOLUTION"])))
    per_energy("angular_resolution", "68 % containment radius (deg)", true, IRF_LIMITS["angular_resolution"],
               lambda d: (*_edges(d["angular"]), d["angular"]["ANGULAR_RESOLUTION_68"].ravel().to_value(u.deg)),
               plot_title="Angular resolution")
    # cuts: bins with at least 10 events (the other ones have the fill value of ctapipe)
    per_energy("theta_cut", "θ cut (deg)", reco, IRF_LIMITS["theta_cut"], lambda d: _cut(d["theta_cuts"]),
               plot_title="Direction cut")
    per_energy("gammaness_cut", "gammaness cut", reco, IRF_LIMITS["gammaness_cut"], lambda d: _cut(d["gh_cuts"]))
    per_energy("background", "rate (s⁻¹ sr⁻¹ TeV⁻¹)", reco, IRF_LIMITS["background"],
               lambda d: (*_edges(d["background"]),
                          d["background"]["BKG"].ravel().to_value(u.Unit("s-1 sr-1 TeV-1"))),
               log=True, plot_title="Background rate")
    per_energy("sensitivity_crab", "flux sensitivity (% Crab)", reco, IRF_LIMITS["crab"],
               lambda d: (*_edges(d["sensitivity"]), 100 * np.ravel(d["sensitivity"]["RELATIVE_SENSITIVITY"])),
               log=True, plot_title="Sensitivity, 50 h", quantity="sensitivity",
               transform=lambda low, high, values: 100 * values / _crab(np.sqrt(low * high)))

    # E² flux sensitivity, with the Crab Nebula
    fig, ax = figure()
    energy = np.geomspace(*IRF_LIMITS["energy"], 50)
    for fraction in (1, 0.1, 0.01):
        ax.plot(energy, fraction * _crab(energy), color=MUTED, lw=1, ls="--")
        ax.annotate(f"{fraction:.0%} Crab", (0.15, fraction * _crab(0.15)), color=MUTED, fontsize=9,
                    xytext=(0, 4), textcoords="offset points")
    for (label, data), color in zip(inputs.items(), COLORS, strict=False):
        row = data["sensitivity"]
        values = row["ENERGY_FLUX_SENSITIVITY"].ravel().to_value(u.TeV / u.cm**2 / u.s)
        _series(ax, *_edges(row), values, color, label if labels else None)
    _reference(ax, reference, "sensitivity")
    ax.set(xscale="log", yscale="log", xlim=IRF_LIMITS["energy"], ylim=IRF_LIMITS["sensitivity"], xlabel=reco,
           ylabel="E² dN/dE (TeV cm⁻² s⁻¹)")
    legend(ax)
    paths.append(save(fig, directory / "sensitivity.png", f"Flux sensitivity, 50 h{suffix}"))

    # energy migration: one figure per input
    for label, data in inputs.items():
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


def plot_gamma_hadron_separation(inputs, directory, prefix="RandomForestClassifier", title="", reference=None):
    """
    Gammaness of the gammas and protons, and ROC curve, one series per input: ``inputs``
    {label: (gamma DL2 file, proton DL2 file)}, ``reference``: table of `read_reference` (ROC)
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
    _reference(ax, reference, "roc")
    ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls="--")
    ax.set(xlim=IRF_LIMITS["fraction"], ylim=IRF_LIMITS["fraction"], xlabel="proton efficiency",
           ylabel="gamma efficiency")
    ax.legend(frameon=False, fontsize=8)
    paths.append(save(fig, directory / "roc.png", f"ROC curve{suffix}"))
    return paths
