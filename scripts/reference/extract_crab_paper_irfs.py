"""
Performance of the SST-1M Crab paper (Alispach et al. 2025, A&A 699, A255, arXiv:2506.01733), extracted
from the vector figures (SVG) of its HTML version into IRF directories of the format of the pipeline
(``ctapipe-compute-irf``: ``irf_point_like.fits.gz``, ``benchmarks_point_like.fits.gz``), for
``sst1mpipe-plot-irfs --reference``: sst1mpipe/resources/reference/crab_paper_2025/<mode>/ze<zenith>/

- tel1 (SST-1M-1 mono), strict_stereo (paper: stereo, events of both telescopes): ze30, ze50,
  benchmarks (energy bias and resolution, 68 % angular resolution, sensitivity)
- tel1, tel2, strict_stereo: ze20, effective area (analysis level)

The quantities the paper does not show (energy dispersion, background, cuts, ...) are not in the files.

    for f in diff_sensitivity energy_res energy_bias angular_res effective_area; do
        curl -sL https://arxiv.org/html/2506.01733v2/$f.svg -o svg/$f.svg
    done
    python scripts/reference/extract_crab_paper_irfs.py svg/ sst1mpipe/resources/reference/crab_paper_2025

The axes are calibrated with the positions of the major ticks of each figure and their values (read on
the figures). Binned figures: one horizontal segment per energy bin, the zenith angle from the style
of the line joining the bins (paper: 30 deg solid, 50 deg dotted).
"""
import re
import sys
from pathlib import Path

import numpy as np

from sst1mpipe.performance.irf_plots import write_reference_directory

SERIES_COLORS = {"#0000ff": "tel1", "#ffa500": "strict_stereo"}
AEFF_COLORS = {"#0000ff": "tel1", "#ffa500": "tel2", "#4daf4a": "strict_stereo"}
AEFF_STAGES = {None: "triggered", "7.4,3.2": "survived cleaning", "12.8,3.2,2,3.2": "quality cuts",
               "2,3.3": "analysis level"}

#: figure: (x major tick values, x log, y major tick values, y log), from the figures
AXES = {
    "diff_sensitivity": ([0.1, 1, 10, 100, 1000], True, [1e-12, 1e-11, 1e-10], True),
    "energy_res": ([0.1, 1, 10, 100], True, list(np.arange(0, 0.85, 0.1)), False),
    "energy_bias": ([0.1, 1, 10, 100], True, list(np.arange(-0.3, 0.35, 0.1)), False),
    "angular_res": ([0.1, 1, 10, 100], True, list(np.arange(0, 0.65, 0.1)), False),
    "effective_area": ([1, 10, 100], True, [1e2, 1e3, 1e4, 1e5, 1e6], True),
}
#: binned figures: quantity
BINNED = {"diff_sensitivity": "sensitivity", "energy_res": "resolution", "energy_bias": "bias",
          "angular_res": "angular_resolution"}
#: energy bins of the binned figures (5 per decade), of the effective area (10 per decade)
ENERGY_BINS = np.logspace(-1, 3, 21)
AEFF_BINS = np.logspace(-1, 3, 41)
REFERENCE = "Alispach et al. 2025, A&A 699, A255, arXiv:2506.01733"


def paths(svg):
    """(stroke color, dash array, list of (x, y)) of the paths of an SVG of matplotlib"""
    out = []
    for attrs in re.findall(r"<path([^>]*)/?>", svg):
        color = re.search(r'stroke="(#[0-9a-f]{6})"', attrs)
        dash = re.search(r'stroke-dasharray="([^"]+)"', attrs)
        d = re.search(r' d="([^"]+)"', attrs)
        if not (color and d):
            continue
        points, x, y = [], None, None
        for command, values in re.findall(r"([MLHVZz])?\s*([-\d.eE ,]*)", d[1]):
            numbers = [float(v) for v in re.findall(r"-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", values)]
            if not numbers:
                continue
            if command == "H":
                x = numbers[-1]
                points.append((x, y))
            elif command == "V":
                y = numbers[-1]
                points.append((x, y))
            else:
                for i in range(0, len(numbers) - 1, 2):
                    x, y = numbers[i], numbers[i + 1]
                    points.append((x, y))
        out.append((color[1], dash[1] if dash else None, points))
    return out


def axis_transforms(svg, figure):
    """Functions path coordinates -> data, from the major ticks (black, 3.5 long)"""
    x_ticks, y_ticks = [], []
    for color, _, points in paths(svg):
        if color != "#000000" or len(points) != 2:
            continue
        (x1, y1), (x2, y2) = points
        if x1 == x2 and abs(abs(y2 - y1) - 3.5) < 1e-3:
            x_ticks.append(x1)
        elif y1 == y2 and abs(abs(x2 - x1) - 3.5) < 1e-3:
            y_ticks.append(y1)
    x_values, x_log, y_values, y_log = AXES[figure]

    def transform(ticks, values, log):
        ticks = sorted(set(round(t, 3) for t in ticks))[: len(values)]
        if len(ticks) != len(values):
            raise ValueError(f"{figure}: {len(ticks)} major ticks for {len(values)} values")
        v = np.log10(values) if log else np.asarray(values)
        slope, intercept = np.polyfit(ticks, v, 1)
        return (lambda p: 10 ** (slope * p + intercept)) if log else (lambda p: slope * p + intercept)

    return transform(x_ticks, x_values, x_log), transform(y_ticks, y_values, y_log)


def binned(svg, figure):
    """Bins (horizontal segments) of each series and zenith angle"""
    to_x, to_y = axis_transforms(svg, figure)
    rows = []
    all_paths = paths(svg)
    for color, regime in SERIES_COLORS.items():
        segments = [p for c, dash, p in all_paths if c == color and dash is None and len(p) == 2
                    and p[0][1] == p[1][1] and p[1][0] > p[0][0]]
        lines = [(dash, p) for c, dash, p in all_paths if c == color and len(p) > 2]
        # contiguous bins: one series
        groups, current = [], [segments[0]]
        for segment in segments[1:]:
            if abs(segment[0][0] - current[-1][1][0]) < 1e-3:
                current.append(segment)
            else:
                groups.append(current)
                current = [segment]
        groups.append(current)
        for group in groups:
            centers = [((s[0][0] + s[1][0]) / 2, s[0][1]) for s in group]
            # the line joining the bins: 30 deg solid, 50 deg dotted
            matches = {}
            for dash, points in lines:
                vertices = np.array(points)
                matches[dash] = sum(np.min(np.hypot(*(vertices - c).T)) < 0.5 for c in centers)
            dash = max(matches, key=matches.get)
            if matches[dash] < len(group) / 2:
                continue  # legend sample
            zenith = 30 if dash is None else 50
            for (x1, y), (x2, _) in group:
                rows.append((regime, zenith, "", to_x(x1), to_x(x2), to_y(y)))
    return rows


def curves(svg, figure, colors, styles):
    """Continuous curves (more than 5 points) of each color and line style"""
    to_x, to_y = axis_transforms(svg, figure)
    rows = []
    for color, dash, points in paths(svg):
        if color in colors and dash in styles and len(points) > 5:
            for x, y in points:
                rows.append((colors[color], styles[dash], to_x(x), to_y(y)))
    return rows


def _per_bin(rows, bins):
    """Values of the bins of the figure in the bins ``bins`` (NaN without value)"""
    values = np.full(len(bins) - 1, np.nan)
    for low, high, value in rows:
        index = np.argmin(np.abs(np.log10(bins[:-1]) - np.log10(low)))
        if np.isclose(np.log10(bins[index]), np.log10(low), atol=0.02) and np.isclose(
                np.log10(bins[index + 1]), np.log10(high), atol=0.02):
            values[index] = value
    return values


def main(svg_directory, output):
    svg_directory, output = Path(svg_directory), Path(output)
    values = {}  # (mode, zenith): {quantity: rows}
    for figure, quantity in BINNED.items():
        for mode, zenith, _, low, high, value in binned((svg_directory / f"{figure}.svg").read_text(), figure):
            values.setdefault((mode, zenith), {}).setdefault(quantity, []).append((low, high, value))
    curves_aeff = {}
    svg = (svg_directory / "effective_area.svg").read_text()
    for mode, stage, x, value in curves(svg, "effective_area", AEFF_COLORS, AEFF_STAGES):
        if stage == "analysis level":
            curves_aeff.setdefault(mode, []).append((x, value))

    comments = {"energy": "true energy; RESOLUTION: (dE/E)_68 of the paper, BIAS: median(E_reco/E_true - 1)",
                "sensitivity": "reconstructed energy; 50 h, 5 sigma",
                "effective_area": "analysis level (gammaness and theta cuts)"}
    for (mode, zenith), quantities in sorted(values.items()):
        write_reference_directory(
            output / mode / f"ze{zenith}", REFERENCE, true_energy_bins=ENERGY_BINS, reco_energy_bins=ENERGY_BINS,
            energy_bias=_per_bin(quantities["bias"], ENERGY_BINS),
            energy_resolution=_per_bin(quantities["resolution"], ENERGY_BINS),
            angular_resolution=_per_bin(quantities["angular_resolution"], ENERGY_BINS),
            sensitivity=_per_bin(quantities["sensitivity"], ENERGY_BINS), comments=comments,
        )
    for mode, points in sorted(curves_aeff.items()):
        x, value = np.array(points).T
        order = np.argsort(x)
        centers = np.sqrt(AEFF_BINS[:-1] * AEFF_BINS[1:])
        inside = (centers >= x.min()) & (centers <= x.max())
        aeff = np.full(len(centers), np.nan)
        aeff[inside] = 10 ** np.interp(np.log10(centers[inside]), np.log10(x[order]), np.log10(value[order]))
        write_reference_directory(output / mode / "ze20", REFERENCE, aeff_energy_bins=AEFF_BINS, effective_area=aeff,
                                  comments=comments)

if __name__ == "__main__":
    main(*sys.argv[1:3])
