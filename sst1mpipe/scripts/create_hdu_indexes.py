#!/usr/bin/env python
"""
Index files of gammapy (``hdu-index.fits.gz``, ``obs-index.fits.gz``) of a selection of DL3 files of
``sst1mpipe-dl2-dl3`` (IRFs in the DL3 files), e.g. several nights of a source for a joint analysis:

    create_hdu_indexes --input-files '/dl3/20240110/*.dl3.fits.gz /dl3/20240111/*.dl3.fits.gz' --output-dir ./crab
"""
import argparse
import glob

from sst1mpipe.io.dl3 import write_index


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Index files of gammapy of DL3 files")
    parser.add_argument(
        "--input-files", "-d", dest="dl3_files", required=True,
        help="DL3 files: glob patterns separated by spaces, e.g. '/dir1/*/*.dl3.fits.gz /dir2/*.dl3.fits.gz'",
    )
    parser.add_argument("--output-dir", "-o", dest="outdir", default="./", help="Directory of the index files")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    files = sorted({path for pattern in args.dl3_files.split() for path in glob.glob(pattern)})
    if not files:
        raise SystemExit(f"No DL3 files: {args.dl3_files}")
    data_store = write_index(files, args.outdir)
    print(data_store.obs_table)


if __name__ == "__main__":
    main()
