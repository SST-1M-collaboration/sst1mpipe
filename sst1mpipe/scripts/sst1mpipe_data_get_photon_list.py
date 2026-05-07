#!/usr/bin/env python

"""
A script to extract a photon list from DL2 applying a simple global 
gammaness cut and event selection for fast theta2 analysis.
- Input is a single DL2 file in HDF format
- Output is hdf file, basicaly DL2 table with gammaness cut and
event selection applied

Usage:

$> python sst1mpipe_data_get_photon_list.py
--input-file SST1M1_20240304_0012_dl2.h5
--output-dir ./
--config sst1mpipe_config.json

"""

import argparse
import logging
import os
import shutil
import sys

from astropy.time import Time

import sst1mpipe
from sst1mpipe.analysis import add_reco_ra_dec
from sst1mpipe.io import check_outdir, load_config, load_dl2_sst1m, write_photon_list
from sst1mpipe.utils import (
    get_horizon_frame,
    get_telescopes,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Extract candidate photon events from data.")

    # Required arguments
    parser.add_argument(
                        '--input-file', '-f',
                        dest='input_file',
                        required=True,
                        help='Path to the DL2 file'
                        )

    parser.add_argument('--config', '-c', action='store', type=str,
                        dest='config_file',
                        help='Path to a configuration file.',
                        required=True
                        )

    # Optional arguments
    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store the output photon lists.',
                        default='./'
                        )
    args = parser.parse_args()
    return args

def main():

    args = parse_args()

    input_file = args.input_file
    outdir = args.outdir
    output_file = os.path.join(outdir, input_file.split('/')[-1].rstrip(".h5") + "_pl.h5")
    output_logfile = os.path.join(outdir, input_file.split('/')[-1].rstrip(".h5") + "_pl.log")
    output_cfgfile = os.path.join(outdir, input_file.split('/')[-1].rstrip(".h5") + "_pl.cfg")

    check_outdir(outdir)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers= [
            logging.FileHandler(output_logfile, 'w+'),
            logging.StreamHandler(stream=sys.stdout)
            ]
    )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    config = load_config(args.config_file)
    shutil.copyfile(args.config_file, output_cfgfile)

    logging.info('Input file: %s', input_file)

    telescopes = get_telescopes(input_file, data_level='dl2')
    tel = telescopes[0]

    dl2_data = load_dl2_sst1m(input_file, tel=tel, config=config, table='pandas', energy_min=0.0)
    
    gammaness_mask = dl2_data.gammaness > config['analysis']['global_gammaness_cut']
    logging.info('Global gammaness cut > %f applied.', config['analysis']['global_gammaness_cut']) 

    dl2_photon_list = dl2_data[gammaness_mask]

    logging.info('N of events of %s after selection cuts: %d', tel, len(dl2_photon_list))

    # transform reco alt az to ra dec
    if tel == 'stereo':
        # For computation of reco ra/dec from reco alt/az in stereo data, it should not matter which telescope frame we use.
        # The respective telescope frames are correctly used in reconstruction itself to get reco alt/az.
        horizon_frame = get_horizon_frame(config=config, telescope='tel_021', times=Time(dl2_photon_list.local_time, format='unix'))
    else:
        horizon_frame = get_horizon_frame(config=config, telescope=tel, times=Time(dl2_photon_list.local_time, format='unix'))

    logging.info('Adding reco RA and DEC ...')
    dl2_photon_list = add_reco_ra_dec(dl2_photon_list, horizon_frame=horizon_frame)

    write_photon_list(
        dl2_photon_list, 
        output_file=output_file, 
        telescope=tel, 
        config=config, 
        )

if __name__ == '__main__':
    main()