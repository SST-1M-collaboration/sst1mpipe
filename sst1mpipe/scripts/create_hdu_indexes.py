#!/usr/bin/env python

"""
A script to create HDU index files for a user-defined selection of DL3 files.
- Inputs are DL3 files to be indexed, together with IRF base directory.
- Outputs are HDU index files

Usage:

$> python create_hdu_indexes.py
--input-files '/tel1/DL3/20240110/*.fits /tel1/DL3/20240111/*.fits'
--irf-dir /IRFs/
--output-dir ./

"""

import argparse
import glob
from sst1mpipe.io import check_outdir
from sst1mpipe.io.dl2_to_dl3 import create_hdu_index
from gammapy.data import DataStore
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Makes index files for all input DL3 files.")

    # Required arguments
    parser.add_argument(
                        '--input-files', '-d',
                        dest='dl3_files',
                        required=True,
                        help='Input DL3 fits files. You may use wild cards, or specify more directories separated with space, e.g. --input-files=\'/dir1/*/CRAB/*fits /dir2/*/CRAB/*fits\''
                        )

    parser.add_argument(
                        '--irf-dir', '-i',
                        dest='irf_dir',
                        default=None,
                        help='Directory with IRFs. This should be just the base directory, specific IRF file in ../data/ is found authomaticaly to all given observing blocks.',
                        )

    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store index files.',
                        default='./'
                        )

    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    dl3_files = args.dl3_files
    irf_dir = args.irf_dir
    outdir = args.outdir

    check_outdir(outdir)

    file_lists = dl3_files.split(' ')
    dl3_file_list = []
    for flist in file_lists:
        dl3_file_list.append(np.array(glob.glob(flist)))
    list_flattened = [item for row in dl3_file_list for item in row]

    #print(list_flattened)
    create_hdu_index(list_flattened, out_dir=outdir, irf_dir=irf_dir)

    data_store = DataStore.from_dir(outdir)
    print(data_store.obs_table)

if __name__ == "__main__":
    main()