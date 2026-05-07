#!/usr/bin/env python

"""
A script to train Random Forests on training MC DL1 files. It can train 
both mono and stereo models. Both are trained per telescope as in the stereo 
reconstruction the resulting parameters are just averaged. For stereo reconstruction, 
however, one may use additional features such as h_max or impact_distance as results 
of geometrical reconstruction.
- Inputs are MC DL1 diffuse gamma file and MC DL1 diffuse proton file. These are usualy 
merged from individual DL1 files using sst1mpipe_merge_hdf5 script to reach satisfactory 
statistics for RF training
- Outputs are trained models in scikit.learn format (.sav)

Usage:

$> python sst1mpipe_mc_train_rfs.py
--input-file-gamma gamma_200_300E3GeV_30_30deg_training_dl1.h5
--input-file-proton proton_400_500E3GeV_30_30deg_training_dl1.h5
--output-dir ./
--config sst1mpipe_config.json
--telescope tel_001
--plot-features
--stereo

"""

import argparse
import logging
import os
import shutil
import sys

import sst1mpipe
from sst1mpipe.io import check_outdir, load_config, load_dl1_sst1m
from sst1mpipe.reco import train_models


def parse_args():
    parser = argparse.ArgumentParser(description="MC train RFs")

    # Required arguments
    parser.add_argument(
                        '--input-file-gamma', '--fg',
                        dest='gammas',
                        required=True,
                        help='Path to the DL1 diffuse gamma file for training'
                        )
    parser.add_argument('--config', '-c', action='store', type=str,
                        dest='config_file',
                        help='Path to a configuration file.',
                        required=True
                        )

    # Optional arguments
    parser.add_argument(
                        '--input-file-proton', '--fp',
                        dest='protons',
                        default=None,
                        help='Path to the DL1 diffuse proton file for training',
                        )
    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store the trained models.',
                        default='./'
                        )

    parser.add_argument(
                        '--telescope', '-t', type=str,
                        dest='telescope',
                        help='Selected telescope: tel_001 or tel_002',
                        default='tel_001'
                        )

    parser.add_argument(
                        '--plot-features',
                        action='store_true',
                        help='Plot feature importances and save figs in the output directory.',
                        dest='plot'
                        )

    parser.add_argument(
                        '--stereo',
                        action='store_true',
                        help='Train RFs for stereo reconstruction.',
                        dest='stereo'
                        )

    args = parser.parse_args()
    return args

def main():

    args = parse_args()

    input_file_gamma = args.gammas
    input_file_proton = args.protons
    outdir = args.outdir
    telescope = args.telescope
    plot = args.plot
    stereo = args.stereo

    check_outdir(outdir)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers= [
            logging.FileHandler(outdir+'/sst1mpipe_rfs_training_'+telescope+'.log', 'w+'),
            logging.StreamHandler(stream=sys.stdout)
            ]
    )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    config = load_config(args.config_file)
    output_cfgfile = os.path.join(outdir, input_file_gamma.split('/')[-1].rstrip(".h5") + "_train.cfg")
    shutil.copyfile(args.config_file, output_cfgfile)

    logging.info('Training Random Forests for Telescope %s', telescope)

    params_gamma = load_dl1_sst1m(input_file_gamma, tel=telescope, config=config, table='pandas', check_finite=True, stereo=stereo, quality_cuts=True)
    if input_file_proton is not None:
        params_protons = load_dl1_sst1m(input_file_proton, tel=telescope, config=config, table='pandas', check_finite=True, stereo=stereo, quality_cuts=True)
    else:
        params_protons = None

    train_models(params_gamma, params_protons, config=config, plot=plot, outdir=outdir, telescope=telescope, stereo=stereo)

if __name__ == '__main__':
    main()
