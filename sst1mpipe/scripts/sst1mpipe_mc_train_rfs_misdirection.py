#!/usr/bin/env python


import sst1mpipe
import argparse
import sys
import os
import logging
import shutil
from sst1mpipe.io import (
    load_dl2_sst1m,
    load_config,
    check_outdir
)
from sst1mpipe.reco import (
    train_rf_misdirection
)

from sst1mpipe.utils import (
    get_telescopes, 
    calculate_misdirection
)

import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="MC train RFs")

    # Required arguments
    parser.add_argument(
                        '--input-file-gamma', '--fg',
                        dest='gammas',
                        required=True,
                        help='Path to the DL2 diffuse gamma file for training. It should be some smaller part of the general training gammas we use in case we do not care about event classes. It should be already reconstructed file, as we need reco direction to calculate misdirection.'
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
                        help='Path to store the trained models.',
                        default='./'
                        )
    parser.add_argument(
                        '--plot-features',
                        action='store_true',
                        help='Plot feature importances and save figs in the output directory.',
                        dest='plot'
                        )

    args = parser.parse_args()
    return args

def main():

    args = parse_args()

    input_file_gamma = args.gammas
    outdir = args.outdir
    plot = args.plot

    check_outdir(outdir)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers= [
            logging.FileHandler(outdir+'/sst1mpipe_rfs_mis_training.log', 'w+'),
            logging.StreamHandler(stream=sys.stdout)
            ]
    )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    config = load_config(args.config_file)
    output_cfgfile = os.path.join(outdir, input_file_gamma.split('/')[-1].rstrip(".h5") + "_train.cfg")
    shutil.copyfile(args.config_file, output_cfgfile)

    tels = get_telescopes(input_file_gamma, data_level="dl2")
    for tel in tels:
        
        dl2_gamma = load_dl2_sst1m(input_file_gamma, tel=tel, config=config, table='pandas')

        # First, we need to calculate misdirection in the DL2 file
        dl2_gamma = calculate_misdirection(dl2_gamma)

        # These features are missing in DL2, we need to add them
        if tel == 'stereo':
            dl2_gamma['log_camera_frame_hillas_intensity_tel1'] = np.log10(dl2_gamma['camera_frame_hillas_intensity_tel1'])
            dl2_gamma['log_camera_frame_hillas_intensity_tel2'] = np.log10(dl2_gamma['camera_frame_hillas_intensity_tel2'])

        logging.info('Training misdirection Random Forests for Telescope %s', tel)
        train_rf_misdirection(dl2_gamma, config=config, outdir=outdir, telescope=tel, plot=plot)
    
if __name__ == '__main__':
    main()