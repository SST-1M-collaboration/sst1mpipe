"""
Performance of the SST-1M Crab paper (Alispach et al. 2025, A&A 699, A255, arXiv:2506.01733), extracted
from the vector figures (SVG) of its HTML version into a table for sst1mpipe-plot-irfs (--reference):
sst1mpipe/resources/reference/sst1m_crab_paper_2025.ecsv.

    for f in diff_sensitivity energy_res energy_bias angular_res effective_area roc; do
        curl -sL https://arxiv.org/html/2506.01733v2/$f.svg -o svg/$f.svg
    done
    python extract_crab_paper_irfs.py svg/ sst1mpipe/resources/reference/sst1m_crab_paper_2025.ecsv

The axes are calibrated with the positions of the major ticks of each figure and their values (read on
the figures). Binned figures: one horizontal segment per energy bin, the zenith angle from the style
of the line joining the bins (paper: 30 deg solid, 50 deg dotted).
"""
import re
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

SERIES_COLORS = {"#0000ff": "SST-1M-1 mono", "#ffa500": "stereo"}
AEFF_COLORS = {"#0000ff": "SST-1M-1 mono", "#ffa500": "SST-1M-2 mono", "#4daf4a": "stereo"}
AEFF_STAGES = {None: "triggered", "7.4,3.2": "survived cleaning", "12.8,3.2,2,3.2": "quality cuts",
               "2,3.3": "analysis level"}

#: figure: (x major tick values, x log, y major tick values, y log), from the figures
AXES = {
    "diff_sensitivity": ([0.1, 1, 10, 100, 1000], True, [1e-12, 1e-11, 1e-10], True),
    "energy_res": ([0.1, 1, 10, 100], True, list(np.arange(0, 0.85, 0.1)), False),
    "energy_bias": ([0.1, 1, 10, 100], True, list(np.arange(-0.3, 0.35, 0.1)), False),
    "angular_res": ([0.1, 1, 10, 100], True, list(np.arange(0, 0.65, 0.1)), False),
    "effective_area": ([1, 10, 100], True, [1e2, 1e3, 1e4, 1e5, 1e6], True),
    "roc": (list(np.arange(0, 1.05, 0.2)), False, list(np.arange(0, 1.05, 0.2)), False),
}
#: binned figures: quantity, unit, x (energy)
BINNED = {
    "diff_sensitivity": ("sensitivity", "TeV cm-2 s-1", "reconstructed energy"),
    "energy_res": ("energy_resolution", "", "true energy"),
    "energy_bias": ("energy_bias", "", "true energy"),
    "angular_res": ("angular_resolution", "deg", "true energy"),
}


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


def curves(svg, figure, colors, styles, min_step=0.0):
    """Continuous curves (more than 5 points) of each color and line style"""
    to_x, to_y = axis_transforms(svg, figure)
    rows = []
    for color, dash, points in paths(svg):
        if color in colors and dash in styles and len(points) > 5:
            last = None
            for i, (x, y) in enumerate(points):
                # points closer than min_step (path units) to the previous one dropped (ROC: ~1000s)
                if last is None or i == len(points) - 1 or np.hypot(x - last[0], y - last[1]) >= min_step:
                    rows.append((colors[color], styles[dash], to_x(x), to_y(y)))
                    last = (x, y)
    return rows


def main(svg_directory, output):
    svg_directory = Path(svg_directory)
    table = Table(names=["quantity", "series", "regime", "zenith", "stage", "x_low", "x_high", "value"],
                  dtype=[str, str, str, int, str, float, float, float])
    for figure, (quantity, _, _) in BINNED.items():
        for regime, zenith, stage, low, high, value in binned((svg_directory / f"{figure}.svg").read_text(), figure):
            table.add_row((quantity, f"{regime}, {zenith} deg", regime, zenith, stage, low, high, value))
    svg = (svg_directory / "effective_area.svg").read_text()
    for regime, stage, x, value in curves(svg, "effective_area", AEFF_COLORS, AEFF_STAGES):
        table.add_row(("effective_area", f"{regime}, 20 deg, {stage}", regime, 20, stage, x, x, value))
    svg = (svg_directory / "roc.svg").read_text()
    for regime, zenith, x, value in curves(svg, "roc", SERIES_COLORS, {None: 30, "1.5,2.475": 50}, min_step=1.5):
        table.add_row(("roc", f"{regime}, {zenith} deg", regime, zenith, "", x, x, value))
    table.meta = {
        "reference": "Alispach et al. 2025, A&A 699, A255, arXiv:2506.01733 (SST-1M Crab Nebula)",
        "extraction": "vector figures of the HTML version (scripts/reference/extract_crab_paper_irfs.py)",
        "units": {"sensitivity": "TeV cm-2 s-1 (E2 dN/dE, 50 h, 5 sigma)", "energy_resolution": "(dE/E)_68",
                  "energy_bias": "median(E_reco / E_true - 1)", "angular_resolution": "deg",
                  "effective_area": "m2", "roc": "x: false positive rate, value: true positive rate"},
        "x": {"sensitivity": "reconstructed energy (TeV)", "energy_resolution": "true energy (TeV)",
              "energy_bias": "true energy (TeV)", "angular_resolution": "true energy (TeV)",
              "effective_area": "true energy (TeV)", "roc": "false positive rate"},
        "figures": {"sensitivity": 11, "energy_resolution": 9, "energy_bias": 9, "angular_resolution": 10,
                    "effective_area": 8, "roc": 7},
    }
    table.write(output, format="ascii.ecsv", overwrite=True)
    for quantity in np.unique(table["quantity"]):
        selected = table[table["quantity"] == quantity]
        print(quantity, {s: int(np.sum(selected["series"] == s)) for s in np.unique(selected["series"])})


if __name__ == "__main__":
    main(*sys.argv[1:3])
