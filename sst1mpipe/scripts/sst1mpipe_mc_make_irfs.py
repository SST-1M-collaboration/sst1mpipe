#!/usr/bin/env python

"""
A script to create IRFs from DL2 MC files.
- Inputs are DL2 point gammas and DL2 diffuse proton files.
- Outputs are IRFs in a format compatible with gammapy.
- Can create biased IRFs in true energy for systematics studies
- Can create separate IRFs for different event classes

Usage:

$> python sst1mpipe_mc_make_irfs.py
--input-file-gamma gamma_200_300E3GeV_30_30deg_testing_dl2.h5
--input-file-proton proton_400_500E3GeV_30_30deg_testing_dl2.h5
--config sst1mpipe_config.json
--output-dir ./
--gammaness-cut-dir ./
--event-classes
--scale-true-energies
"""

import argparse
import sys

from sst1mpipe.utils import get_telescopes
from sst1mpipe.io import check_outdir, load_dl2_sst1m, load_config
from sst1mpipe.performance import irf_maker
import sst1mpipe
import logging
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="IRF maker")

    # Required arguments
    parser.add_argument(
                        '--input-file-gamma', '-fg',
                        dest='gammas',
                        required=True,
                        help='Path to the DL2 testing diffuse gamma file'
                        )

    parser.add_argument(
                        '--input-file-proton', '-fp',
                        dest='protons',
                        default=None,
                        help='Path to the DL2 testing diffuse proton file.',
                        )

    parser.add_argument('--config', '-c', action='store', type=str,
                        dest='config_file',
                        help='Path to a configuration file.',
                        required=True
                        )

    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store the IRFs. This should be just base directory, the rest is created authomaticaly.',
                        default='./IRFs/'
                        )
    
    parser.add_argument(
                        '--gammaness-cut-dir', type=str,
                        dest='gammaness_cut_dir',
                        help='Directory with gammaness cut tables optimized on MC. If empty, global gammaness cut from the config file is used.',
                        default=None
                        )
    parser.add_argument(
                        '--scale-true-energies',
                        action='store_true',
                        help='Enables scaling of true energy for production of a biased IRFs. If used, a global scale on true energy of each shower is applied from config[\"analysis\"][\"true_energy_scaling_factor\"]. Note that this does not modify the true energy values stored in the DL2 files.',
                        dest='true_energy_scaling'
                        )
    parser.add_argument(
                        '--event-classes',
                        action='store_true',
                        dest='event_classes',
                        help='If used, IRFs are created for each event class found in the input DL2 files.',
                        )           
    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    input_file_gamma = args.gammas
    input_file_proton = args.protons
    outdir = args.outdir
    gammaness_cut_dir = args.gammaness_cut_dir
    true_energy_scaling = args.true_energy_scaling
    event_classes = args.event_classes

    check_outdir(outdir)

    # This is not very nice, but I am not yet sure how to put the log file in the final directory, 
    # containing only outputs for given telescope, if we want to iterate over them below
    zenith = input_file_gamma.split('/')[-1].split('_')[3]
    logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers= [
                logging.FileHandler(outdir+'/sst1mpipe_irfs_'+str(sst1mpipe.__version__)+'_' + zenith + '.log', 'w+'),
                logging.StreamHandler(stream=sys.stdout)
                ]
        )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    telescopes = get_telescopes(input_file_gamma, data_level="dl2")

    for tel in telescopes:

        if gammaness_cut_dir is not None:
            gammaness_cuts = glob.glob(gammaness_cut_dir + '/gammaness_cuts*' + tel + '*h5')[0]
            logging.info('Gammaness cut file used: {}'.format(gammaness_cuts))
        else:
            gammaness_cuts = None

        if event_classes:
            config = load_config(args.config_file)
            dl2_mc_gamma  = load_dl2_sst1m(input_file_gamma,
                                       config=config,
                                       tel=tel,
                                       table='pandas')
            evttypes = dl2_mc_gamma['event_type'].unique()

        else: evttypes = []

        if len(evttypes):
            for event_class in evttypes:

                maker = irf_maker(config_filename = args.config_file,
                                    mc_gamma_filename  = input_file_gamma,
                                    mc_proton_filename = input_file_proton,
                                    mc_tel_setup = tel,
                                    point_like_offset = None,
                                    output_dir = outdir,
                                    gammaness_cuts = gammaness_cuts,
                                    true_energy_scaling = true_energy_scaling,
                                    event_class = event_class,
                                    )

                maker.make_all_irfs()
        else:
                maker = irf_maker(config_filename = args.config_file,
                                    mc_gamma_filename  = input_file_gamma,
                                    mc_proton_filename = input_file_proton,
                                    mc_tel_setup = tel,
                                    point_like_offset = None,
                                    output_dir = outdir,
                                    gammaness_cuts = gammaness_cuts,
                                    true_energy_scaling = true_energy_scaling,
                                    )

                maker.make_all_irfs()

if __name__ == "__main__":
    main()


       
