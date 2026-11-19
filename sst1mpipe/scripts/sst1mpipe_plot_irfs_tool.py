"""
``sst1mpipe-plot-irfs``: plots of the IRFs and benchmarks of ``ctapipe-compute-irf`` (directories of
``scripts/camk/mc_irf.bash``) and of the gamma/hadron separation of the DL2, one plot per figure with
the same axis limits (`sst1mpipe.ml.irf_plots`). With several inputs, one series per input (e.g. the
modes tel1, tel2, stereo, strict_stereo).

    sst1mpipe-plot-irfs --irf tel1=irf/optimized/tel1/ze20 --dl2 tel1=dl2/ze20/tel1 -o plots/tel1
    sst1mpipe-plot-irfs --irf tel1=irf/optimized/tel1/ze20 --irf stereo=irf/optimized/stereo/ze20 -o plots/comparison

With ``--reference``, the performance of the SST-1M Crab paper (arXiv:2506.01733: sensitivity, energy
and angular resolution, energy bias, effective area, ROC; or another file of the same format) is drawn
with the series, e.g. at 30 deg (paper: 30 and 50 deg; effective area: 20 deg):

    sst1mpipe-plot-irfs --irf stereo=irf/paper/strict_stereo/ze20 --reference --reference-zenith 30 -o plots/paper
"""
import argparse
from pathlib import Path

import numpy as np

from sst1mpipe.ml.irf_plots import (
    REFERENCE_FILE,
    plot_gamma_hadron_separation,
    plot_irfs,
    read_irf_directory,
    read_reference,
)


def labelled(text):
    label, _, path = text.partition("=")
    if not path:
        raise argparse.ArgumentTypeError(f"expected <label>=<directory>, got {text}")
    return label, Path(path)


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--irf", type=labelled, action="append", required=True,
                        help="<label>=<directory with irf_point_like.fits.gz, benchmarks_point_like.fits.gz, cuts.fits>")
    parser.add_argument("--dl2", type=labelled, action="append", default=[],
                        help="<label>=<directory with gamma_point.dl2.h5 and proton_test.dl2.h5> (gammaness, ROC)")
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument("--title", default="", help="added to the titles of the figures")
    parser.add_argument("--reference", nargs="?", const=REFERENCE_FILE, type=Path, default=None,
                        help=f"reference performance drawn with the series (default file: {REFERENCE_FILE.name},"
                             " the SST-1M Crab paper)")
    parser.add_argument("--reference-zenith", type=int, action="append", default=[],
                        help="only the reference series of this zenith angle (deg; repeatable; default: all)."
                             " The effective area of the paper (20 deg) is always drawn")
    parser.add_argument("--reference-series", action="append", default=[],
                        help="only the reference series containing this text (repeatable), e.g. stereo")
    options = parser.parse_args(args)

    reference = None
    if options.reference is not None:
        reference = read_reference(options.reference)
        keep = reference["quantity"] == "effective_area"
        if options.reference_zenith:
            keep |= np.isin(reference["zenith"], options.reference_zenith)
        else:
            keep[:] = True
        if options.reference_series:
            keep &= [any(text in series for text in options.reference_series) for series in reference["series"]]
        reference = reference[keep]
    paths = plot_irfs({label: read_irf_directory(path) for label, path in options.irf}, options.output_dir,
                      title=options.title, reference=reference)
    if options.dl2:
        inputs = {label: (path / "gamma_point.dl2.h5", path / "proton_test.dl2.h5") for label, path in options.dl2}
        paths += plot_gamma_hadron_separation(inputs, options.output_dir, title=options.title, reference=reference)
    print(f"{len(paths)} figures in {options.output_dir}")


if __name__ == "__main__":
    main()
