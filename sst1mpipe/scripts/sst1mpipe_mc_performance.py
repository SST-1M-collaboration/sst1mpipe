#!/usr/bin/env python

"""
A script to get SST-1M performance from DL2 MC files.
- Inputs are DL2 point gammas and DL2 diffuse proton files.
- Outputs are hdf tables and plots for energy and angular resolution, sensitivity and ROC curve. 
It can also print a time needed for detection of a source of selection. If application of energy 
dependent gammaness cut is selected, the cuts are stored in separate hdf file. This can be used
to produce IRFs and DL3 files.

Usage:

$> python sst1mpipe_mc_performance.py
--input-file-gamma gamma_200_300E3GeV_30_30deg_testing_dl2.h5
--input-file-proton proton_400_500E3GeV_30_30deg_testing_dl2.h5
--config sst1mpipe_config.json
--output-dir ./
--gammaness-cuts efficiency
--source-detection CRAB_MAGIC_JHEAP2015
--energy-min 0.0
--save-fig
--save-hdf
--rf-performance
--sensitivity
--gamma-off-correction

"""

import sst1mpipe
import os
import sys
import argparse
import numpy as np
import astropy.units as u
import shutil
import logging
from sst1mpipe.utils import (
    get_telescopes
)
from sst1mpipe.io import (
    load_config,
    check_outdir
)
from sst1mpipe.performance import (
    evaluate_performance,
    sensitivity
)

def parse_args():
    parser = argparse.ArgumentParser(description="Evaluation of performance from reconstructed MC DL2 files")

    # Required arguments
    parser.add_argument(
                        '--input-file-gamma', '-fg',
                        dest='gammas',
                        required=True,
                        help='Path to the DL2 testing point gamma file'
                        )
    parser.add_argument('--config', '-c', action='store', type=str,
                        dest='config_file',
                        help='Path to a configuration file.',
                        required=True
                        )

    # Optional arguments
    parser.add_argument(
                        '--input-file-proton', '-fp',
                        dest='protons',
                        default=None,
                        help='Path to the DL2 testing diffuse proton file. Not necessary if one wants to estimate only angular or energy resolution.',
                        )
    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store the output figures',
                        default='./'
                        )
    parser.add_argument(
                        '--save-fig',
                        action='store_true',
                        help='Save figures in some default ctaplot formating',
                        dest='save_fig'
                        )
    parser.add_argument(
                        '--save-hdf',
                        action='store_true',
                        help='Save performance in hdf tables.',
                        dest='save_hdf'
                        )
    parser.add_argument(
                        '--gammaness-cuts',
                        action='store', type=str,
                        help='Select method of gammaness cuts to be applied: \'global\' (default), \'efficiency\', \'significance\'.',
                        dest='gammaness_cuts',
                        default='global'
                        )
    parser.add_argument(
			'--theta2-cuts',
			action='store', type=str,
			help='Select method of theta2 cuts to be applied: \'global\' (default), \'efficiency\'.',
			dest='theta2_cuts',
			default='global'
			)
    parser.add_argument(
                        '--rf-performance',
                        action='store_true',
                        help='Calculate energy and angular resolution, and ROC curves.',
                        dest='rf_performance',
                        default=False
                        )
    parser.add_argument(
                        '--sensitivity',
                        action='store_true',
                        help='Calculate sensitivity.',
                        dest='sensitivity',
                        default=False
                        )
    parser.add_argument(
                        '--source-detection', '--s',
                        dest='source_detection',
                        default='',
                        help='Source of interest. Time to the source detection is calculated. Spectrum must be specified in performance.spectra, e.g. CRAB_HEGRA, CRAB_MAGIC_JHEAP2015...',
                        )
    parser.add_argument(
                        '--energy-min',
                        action='store', type=float,
                        help='Additional cut on minimum reconstructed energy in TeV. Sensitivity is not very reliable around the energy threshold.',
                        dest='energy_min',
                        default=0.0
                        )
    parser.add_argument(
                        '--gamma-off-correction',
                        action='store_true',
                        help='Takes into account a number of badly reconstructed gammas which fall in the OFF region for the sensitivity calculations. These contribute to the background rate. This option so far only works for MONO.',
                        dest='gamma_off_correction',
                        default=False
                        )
    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    input_file_gamma = args.gammas
    input_file_proton = args.protons
    outdir = args.outdir
    save_fig = args.save_fig
    save_hdf = args.save_hdf
    gammaness_cuts = args.gammaness_cuts
    theta2_cuts = args.theta2_cuts
    get_rf_performance = args.rf_performance
    get_sensitivity = args.sensitivity
    source_detection = args.source_detection
    energy_min = args.energy_min
    gamma_off = args.gamma_off_correction

    if len(source_detection) > 0 and not get_sensitivity:
        logging.error('Time to source detection is calculated in the sensitivity module. It cannot be calculated alone. Use the switch \'--sensitivity\' as well.')
        exit()

    check_outdir(outdir)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers= [
            logging.FileHandler(outdir+'/sst1mpipe_performance.log', 'w+'),
            logging.StreamHandler(stream=sys.stdout)
            ]
    )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    config = load_config(args.config_file)
    telescopes = get_telescopes(input_file_gamma, data_level="dl2")

    # We want to save config file with cuts in the output directory with performance
    output_cfgfile = os.path.join(outdir, input_file_gamma.split('/')[-1].rstrip(".h5") + "_analysis.cfg")
    shutil.copyfile(args.config_file, output_cfgfile)

    logging.info('Input file gamma: %s', input_file_gamma)
    if input_file_proton is not None:
        logging.info('Input file proton: %s', input_file_proton)

    for tel in telescopes:

        if get_rf_performance:
            evaluate_performance(gamma_file=input_file_gamma, proton_file=input_file_proton, outdir=outdir, config=config, telescope=tel, save_fig=save_fig, save_hdf=save_hdf, gammaness_cuts=gammaness_cuts)
        if get_sensitivity:
            if input_file_proton is not None:
                sensitivity(input_file_gamma, input_file_proton, outdir=outdir, config=config, telescope=tel, save_fig=save_fig, save_hdf=save_hdf, gammaness_cuts=gammaness_cuts, theta2_cuts=theta2_cuts, source_detection=source_detection, energy_min=energy_min, gamma_off=gamma_off)
            else:
                logging.error('DL2 proton file not specified. Sensitivity cannot be estimated.')
        if tel == "stereo":
            break

if __name__ == '__main__':
    main()
