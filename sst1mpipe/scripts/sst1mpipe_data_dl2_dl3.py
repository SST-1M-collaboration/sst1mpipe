#!/usr/bin/env python

"""
A script to create DL3 data files from DL2 data files. It merges DL2 HDF 
per-run files into per-wobble DL3 fits files. It also finds proper IRF
based on zenith, azimuth and NSB level. It also creates per-night index files.
- Input is a directory with DL2 files and directory with IRFs
- Outputs are DL3 files and hdu indexes
- If list of event classes is provided, the DL3 files are calculated for all of them and stored in separate subdirs event_class_X

Usage:

$> python sst1mpipe_data_dl2_dl3.py
--input-dir /data/dl2/
--output-dir ./
--config sst1mpipe_config.json
--irf-dir /IRFs/
--target-name Crab
--target-ra 85.0
--target-dec 25.0
--gammaness-cut-dir /gammaness_cuts/
--event-classes

"""

import numpy as np
import argparse
import sys
import os
import logging

import sst1mpipe
from sst1mpipe.io import check_outdir
from sst1mpipe.utils import get_target_pos

from sst1mpipe.io.dl2_to_dl3 import (
    dl2_dir_to_dl3,
    create_hdu_index
    )

from gammapy.data import DataStore


def parse_args():
    parser = argparse.ArgumentParser(description="Makes DL3 from all DL2 files in given directory.")

    parser.add_argument(
                        '--input-dir', '-d',
                        dest='dl2_dir',
                        required=True,
                        help='Directory with input DL2 files.'
                        )

    parser.add_argument(
                        '--irf-dir', '-i',
                        dest='irf_dir',
                        default=None,
                        help='Directory with IRFs. This should be just the base directory, specific IRF file in ../data/ is found authomaticaly to all given observing blocks.',
                        required=True
                        )

    parser.add_argument('--config', '-c', action='store', type=str,
                        dest='config_file',
                        help='Path to a configuration file.',
                        required=True
                        )

    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store DL3 files.',
                        default='./'
                        )

    parser.add_argument(
                        '--target-name', type=str,
                        dest='target_name',
                        help='Name of the target, e.g. Crab',
                        default=None
                        )

    parser.add_argument(
                    '--target-ra', type=float,
                    dest='target_ra',
                    help='Target RA (deg)',
                    default=None
                    )

    parser.add_argument(
                    '--target-dec', type=float,
                    dest='target_dec',
                    help='Target DEC (deg)',
                    default=None
                    )

    parser.add_argument(
                    '--gammaness-cut-dir', type=str,
                    dest='gammaness_cuts_dir',
                    help='Dir with gammaness cuts optimized on MC. The subdir structure should follow the same logic as the RF model directories. By default global gammaness cut from the config file is used.',
                    default=None
                    )
    parser.add_argument(
                    '--event-classes', type=int,
                    dest='event_classes',
                    nargs='+',
                    default=[],
                    help='Event classes to process, e.g. --event-classes 1 2 3 4. If used, DL3 files are created for each event class in the list (if proper IRFs are found). If not used, event classes are ignored completely.')

    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    dl2_dir = args.dl2_dir
    irf_dir = args.irf_dir
    outdir = args.outdir
    target_name = args.target_name
    target_ra = args.target_ra
    target_dec = args.target_dec
    gammaness_cuts_dir = args.gammaness_cuts_dir
    event_classes = args.event_classes

    check_outdir(outdir)

    logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers= [
                logging.FileHandler(outdir+'/sst1mpipe_dl2_dl3_'+str(sst1mpipe.__version__) + '.log', 'w+'),
                logging.StreamHandler(stream=sys.stdout)
                ]
        )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    target_pos, target_name = get_target_pos(target_name=target_name, ra=target_ra, dec=target_dec)

    dl3_files = dl2_dir_to_dl3(
                target_name = target_name,
                target_pos = target_pos,
                dl2_dir       = dl2_dir,
                config_file   = args.config_file,
                out_dir       = outdir,
                gammaness_cuts = gammaness_cuts_dir,
                event_classes=event_classes
                )
    
    if len(event_classes):
        for evtclass in event_classes:
            outdir_c = outdir + '/event_class_' + str(evtclass)
            create_hdu_index(dl3_files['event_class_'+str(evtclass)], out_dir=outdir_c, irf_dir=irf_dir)
            data_store = DataStore.from_dir(outdir_c)
            print(data_store.obs_table)
    else:
        create_hdu_index(dl3_files['all'], out_dir=outdir, irf_dir=irf_dir)
        data_store = DataStore.from_dir(outdir)
        print(data_store.obs_table)

if __name__ == "__main__":
    main()