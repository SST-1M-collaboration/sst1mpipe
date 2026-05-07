#!/usr/bin/env python

"""
A script to merge MC DL1 files in order to prepare a single (per particle) file 
for RF training/testing.
- Input is a directory containing DL1 MC files.
- Output is a merged DL1 file containing all events from the input files.

Usage:

$> python sst1mpipe_merge_hdf5.py
--input-dir /mc_prod/diffuse_gammas/dl1/
--output-file merged_dl1_diffuse_gammas.h5
--pattern gammas*.h5
--skip-images
--overwrite
--skip-broken

"""

import argparse

from ctapipe.tools.merge import MergeTool


def parse_args():
    parser = argparse.ArgumentParser(description="MC/Data merge hdf5 tables")

    # Required arguments
    parser.add_argument(
                        '-d', '--input-dir',
                        help='path to the source directory of files',
                        required=True,
                        dest='input_dir'
                        )

    # Optional arguments
    parser.add_argument(
                        '-o', '--output-file',
                        help='Path of the resulting merged file',
                        default='./merge.h5',
                        dest='output_file'
                    )

    parser.add_argument(
                        '--skip-images',
                        action='store_true',
                        help='Do not include images in output file',
                        dest='no_images'
                        )

    parser.add_argument(
                        '--overwrite',
                        action='store_true',
                        help='Overwrite the output file if already exists',
                        dest='overwrite'
                        )

    parser.add_argument(
                        '--pattern',
                        help='Specific file pattern for the input files',
                        default='*.h5',
                        dest='pattern'
                    )

    parser.add_argument(
                        '--skip-broken',
                        action='store_true',
                        help='Skip files that cannot be merged instead of raising an error',
                        dest='skip_broken'
                    )

    args = parser.parse_args()
    return args

def main():

    args = parse_args()

    input_dir = args.input_dir
    output_file = args.output_file
    no_images = args.no_images
    overwrite = args.overwrite
    pattern = args.pattern
    skip_broken_files = args.skip_broken

    merge = MergeTool(input_dir=input_dir, output_path=output_file, skip_broken_files=skip_broken_files, skip_images=no_images, skip_simu_images=no_images, overwrite=overwrite, progress_bar=True, file_pattern=pattern)
    merge.run()

if __name__ == '__main__':
    main()
