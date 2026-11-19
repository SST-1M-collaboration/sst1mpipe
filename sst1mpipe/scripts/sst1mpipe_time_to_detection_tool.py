"""
``sst1mpipe-time-to-detection``: time to detect point sources at 5 sigma with the point-like IRFs of
``ctapipe-compute-irf`` (`sst1mpipe.performance.detection`), above minimum reconstructed energies, with
the conditions of the sensitivity of ctapipe (at least 10 signal events, signal at least 5 % of the
background), 5 OFF regions.

    sst1mpipe-time-to-detection --irf irf/optimized/stereo/ze20 --source CRAB_HEGRA --source MRK421_VERITAS_LOW \\
        --energy-min 1 --energy-min 3 -o time_to_detection.png

Sources: the spectra of `sst1mpipe.performance.spectra` (``--list-sources``). The minimum energies are
rounded to the lower edges of the reconstructed energy bins of the background IRF. With ``-o``, the figure
of the time to detection versus the minimum energy.
"""
import argparse
from pathlib import Path

import astropy.units as u
import numpy as np

from sst1mpipe.performance import spectra
from sst1mpipe.performance.detection import (
    detection_rates,
    source_spectrum,
    time_to_detection,
    time_to_detection_table,
)


def irf_file(path):
    path = Path(path)
    return path / "irf_point_like.fits.gz" if path.is_dir() else path


def plot(irf, sources, output, n_off, title):
    from sst1mpipe.ml.plots import COLORS, MUTED, figure, save
    from sst1mpipe.performance.irf_plots import IRF_LIMITS

    fig, ax = figure()
    for source, color in zip(sources, COLORS * 2, strict=False):
        rates = detection_rates(irf, source_spectrum(source))
        energies = rates.energy_edges[:-1]
        times = u.Quantity([time_to_detection(*rates.above(e), n_off=n_off) for e in energies]).to_value(u.h)
        ax.step(energies.to_value(u.TeV), np.where(np.isfinite(times), times, np.nan), where="post", color=color, lw=2,
                label=source)
    for hours in (1, 10, 50):
        ax.axhline(hours, color=MUTED, lw=1, ls=":")
    ax.set(xscale="log", yscale="log", xlim=IRF_LIMITS["energy"], ylim=(0.01, 1000),
           xlabel="minimum reconstructed energy (TeV)", ylabel="time to 5σ (h)")
    ax.legend(frameon=False, fontsize=8)
    save(fig, output, title or "Time to detection")


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--irf", help="irf_point_like.fits.gz of ctapipe-compute-irf, or its directory")
    parser.add_argument("--source", action="append", default=[],
                        help="spectrum of sst1mpipe.performance.spectra (repeatable, default CRAB_HEGRA)")
    parser.add_argument("--energy-min", action="append", type=float, default=[],
                        help="minimum reconstructed energy in TeV (repeatable, default 0.3, 1, 3, 10)")
    parser.add_argument("--n-off", type=int, default=5, help="number of OFF regions")
    parser.add_argument("-o", "--output", type=Path, help="figure of the time to detection versus the minimum energy")
    parser.add_argument("--title", default="", help="title of the figure")
    parser.add_argument("--list-sources", action="store_true", help="list the sources and exit")
    options = parser.parse_args(args)
    if options.list_sources:
        print("\n".join(name for name in spectra.__all__ if "UNIT" not in name and name != "DAMPE_P_He_SPECTRUM"))
        return
    if options.irf is None:
        parser.error("--irf is required")
    irf = irf_file(options.irf)
    sources = options.source or ["CRAB_HEGRA"]
    table = time_to_detection_table(irf, sources, options.energy_min or [0.3, 1, 3, 10], n_off=options.n_off)
    for name in ("signal_rate", "background_rate", "time", "energy_min"):
        table[name].info.format = ".3g"
    table.pprint(max_lines=-1, max_width=-1)
    if options.output:
        plot(irf, sources, options.output, options.n_off, options.title)
        print(f"Figure: {options.output}")


if __name__ == "__main__":
    main()
