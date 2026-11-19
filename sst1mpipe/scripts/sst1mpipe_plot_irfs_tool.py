"""
``sst1mpipe-plot-irfs``: plots of the IRFs, benchmarks and cuts of ``ctapipe-compute-irf`` and
``ctapipe-optimize-event-selection`` and of the gamma/hadron separation of the DL2, one plot per figure
with the same axis limits (`sst1mpipe.performance.irf_plots`). With several inputs, one series per input (e.g. the
modes tel1, tel2, stereo, strict_stereo).

    sst1mpipe-plot-irfs --irf tel1=irf/optimized/tel1/ze20 --dl2 tel1=dl2/ze20/tel1 -o plots/tel1
    sst1mpipe-plot-irfs --irf tel1=irf/optimized/tel1/ze20 --irf stereo=irf/optimized/stereo/ze20 -o plots/comparison

``--reference``: other IRF directories of the same format, drawn as lines, e.g. the published
performance with ``ref:<directory>`` (in `sst1mpipe.performance.irf_plots.REFERENCE_DIRECTORY`):
the SST-1M Crab paper (arXiv:2506.01733), ``crab_paper_2025/<mode>/ze<zenith>``: tel1, strict_stereo
at ze30, ze50 (energy bias and resolution, angular resolution, sensitivity), tel1, tel2, strict_stereo
at ze20 (effective area); the sensitivity of ``lhaaso_1yr``, ``hawc_1yr``, ``veritas_50h``,
``hess_50h``, ``astri_50h``

    sst1mpipe-plot-irfs --irf strict_stereo=irf/paper/strict_stereo/ze20 \\
        --reference "paper 30°=ref:crab_paper_2025/strict_stereo/ze30" \\
        --reference "paper 20°=ref:crab_paper_2025/strict_stereo/ze20" --reference "ASTRI=ref:astri_50h" -o plots/paper
"""
import argparse
from pathlib import Path

from sst1mpipe.performance.irf_plots import REFERENCE_DIRECTORY, plot_gamma_hadron_separation, plot_irfs, read_irf_directory


def labelled(text):
    label, _, path = text.partition("=")
    if not path:
        raise argparse.ArgumentTypeError(f"expected <label>=<directory>, got {text}")
    if path.startswith("ref:"):
        return label, REFERENCE_DIRECTORY / path.removeprefix("ref:")
    return label, Path(path)


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--irf", type=labelled, action="append", default=[],
                        help="<label>=<directory with irf_point_like.fits.gz, benchmarks_point_like.fits.gz, cuts.fits>")
    parser.add_argument("--dl2", type=labelled, action="append", default=[],
                        help="<label>=<directory with gamma_point.dl2.h5 and proton_test.dl2.h5> (gammaness, ROC)")
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument("--title", default="", help="added to the titles of the figures")
    parser.add_argument("--reference", type=labelled, action="append", default=[],
                        help="<label>=<IRF directory>, drawn as lines; ref:<directory>: published performance, e.g. ref:crab_paper_2025/tel1/ze30")
    options = parser.parse_args(args)
    if not options.irf and not options.reference:
        parser.error("no --irf nor --reference")

    paths = plot_irfs({label: read_irf_directory(path) for label, path in options.irf}, options.output_dir,
                      title=options.title,
                      references={label: read_irf_directory(path) for label, path in options.reference})
    if options.dl2:
        inputs = {label: (path / "gamma_point.dl2.h5", path / "proton_test.dl2.h5") for label, path in options.dl2}
        paths += plot_gamma_hadron_separation(inputs, options.output_dir, title=options.title)
    print(f"{len(paths)} figures in {options.output_dir}")


if __name__ == "__main__":
    main()
