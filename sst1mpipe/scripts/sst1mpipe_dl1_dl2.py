#!/usr/bin/env python

"""
A script to reconstruct events with DL1 information using pre-trained RF models. 
Works on both MC and data, mono and stereo.
- Inputs is a single DL1 file in HDF format (output of sst1mpipe_r0_dl1 or 
sst1mpipe_dl1_dl1_stereo)
- Output is hdf file with a table of DL2 parameters (energy, gammaness, src_x, 
src_y, ra, dec). This file still contains DL1 table for convenience.

Usage:

$> python sst1mpipe_dl1_dl2.py
--input-file SST1M1_20240304_0012_dl1.h5
--output-dir ./
--config sst1mpipe_config.json
--models-dir ./rf_models/
--stereo

"""

import sst1mpipe
import os
import sys
import argparse
import shutil
import pandas as pd
import logging
from sst1mpipe.reco import (
    apply_models,
    stereo_reconstruction
)
from sst1mpipe.utils import (
    get_telescopes,
    energy_min_cut,
    check_mc,
    get_closest_rf_model
)
from sst1mpipe.io import (
    load_dl1_sst1m, 
    write_dl2,
    load_config,
    check_outdir,
    write_dl2_info
)

def parse_args():
    parser = argparse.ArgumentParser(description="MC/Data DL1 to DL2")

    # Required arguments
    parser.add_argument(
                        '--input-file', '-f',
                        dest='input_file',
                        required=True,
                        help='Path to the DL1 file to be reconstructed'
                        )

    parser.add_argument('--config', '-c', action='store', type=str,
                        dest='config_file',
                        help='Path to a configuration file.',
                        required=True
                        )

    parser.add_argument(
                        '--models-dir', '-m', type=str,
                        dest='models_dir',
                        help='Path to stored trained RFs, or path to general production directory where subdirectories with models in the following naming format are expected: \'zeXX_azXX_nsbXX\'',
                        required=True
                        )

    # Optional arguments
    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store the output DL2 file',
                        default='./'
                        )

    parser.add_argument(
                        '--stereo',
                        action='store_true',
                        help='Stereo reconstruction. It assumes that there are coincident events in DL1 identified by \'event id\', meaning that each \'event id\' is exactly two times in the DL1 file.',
                        dest='stereo'
                        )
    parser.add_argument(
                        '--scale-intensities',
                        action='store_true',
                        help='Enables Hillas intensity scaling on the DL1 level. If used, for the reconstrucion, a global scale on Hillas Intensity of each shower is applied from config[\"NsbCalibrator\"][\"intensity_correction\"]. Note that this does not modify the Hillas intensity stored in the DL1 table of the DL2 file.',
                        dest='scale_int'
                        )

    args = parser.parse_args()
    return args

def main():

    args = parse_args()

    input_file = args.input_file
    outdir = args.outdir
    output_file = os.path.join(outdir, input_file.split('/')[-1].rstrip(".h5") + "_dl2.h5")
    output_logfile = os.path.join(outdir, input_file.split('/')[-1].rstrip(".h5") + "_dl2.log")
    output_cfgfile = os.path.join(outdir, input_file.split('/')[-1].rstrip(".h5") + "_dl2.cfg")
    models_dir = args.models_dir
    stereo = args.stereo
    scale_int = args.scale_int

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

    telescopes = get_telescopes(input_file)

    logging.info('Input file: %s', input_file)

    # DL2 file is just a copy of DL1 file with DL2 table added. We could consider removal of DL1 table to save some disk space, 
    # but in that case, we would loose the information about all triggered events
    # because we cannot keep NaNs in the DL2 table as RFs are not able to reconstruct these and the reconstruction fails
    shutil.copyfile(input_file, output_file)
    
    ismc = check_mc(input_file)

    # Cut on minimum mc_energy in the output file, which is needed if we want to safely combine MC from different productions
    # NOTE: This doesn't change the mc and histogram tab in the output files and this must be taken care of in performance
    # evaluation. In ideal world, this shouldn't have any effect, because DL1 would already be produced with this cut
    # applied.
    if ismc:
        energy_min_cut(output_file, config=config)

    if stereo:
        dl2_0 = pd.DataFrame()
        for tel in telescopes:

            logging.info('Reconstruction for %s', tel)
            dl1 = load_dl1_sst1m(output_file, tel=tel, config=config, table='pandas', stereo=stereo, check_finite=True, quality_cuts=False, scale_intensities=scale_int)
            # JJ: I turned off selection based on NSB as it is probably not going to be used anyway
            """
            if 'meanQ' in dl1:
                meanQ = np.average(dl1['meanQ'])
                models_dir_auto = get_closest_rf_model(dl1, models_dir=models_dir, nsb_level=meanQ, tel=tel, config=config)
            else:
            """
            models_dir_auto = get_closest_rf_model(dl1, models_dir=models_dir)
            try:
                dl2_0 = pd.concat([dl2_0, apply_models(dl1, models_dir=models_dir_auto, config=config, telescope=tel, stereo=stereo, mc=ismc)])
            except:
                logging.error('RF application failed.')
                os.remove(output_file)
                exit()

        # We have reconstructed energies, gammaness and disp norm per telescope.
        # Now we need to group everything per event and calculate averages.
        # For arrival dirrection, we apply MARS like reconstruction
        #try:
        dl2 = stereo_reconstruction(dl2_0, config=config, ismc=ismc, telescopes=telescopes)
        #except:
        #    logging.error('Stereo reconstruction failed.')
        #    os.remove(output_file)
        #    exit()

        # NOTE: DL2 table from stereo is stored in dl2/event/parameters/stereo
        # This is probably not in accordance with the ctapipe datamodel and should be changed in the future
        try:
            write_dl2(dl2, output_file=output_file, telescope='stereo', config=config, mode='a')
        except:
            logging.error('Writting DL2 file failed.')
            os.remove(output_file)
            exit()

    else:
        for tel in telescopes:

            logging.info('Reconstruction for %s', tel)
            dl1 = load_dl1_sst1m(output_file, tel=tel, config=config, table='pandas', stereo=stereo, check_finite=True, quality_cuts=False, scale_intensities=scale_int)
            # JJ: I turned off selection based on NSB as it is probably not going to be used anyway
            """
            if 'meanQ' in dl1:
                meanQ = np.average(dl1['meanQ'])
                models_dir_auto = get_closest_rf_model(dl1, models_dir=models_dir, nsb_level=meanQ, tel=tel, config=config)
            else:
            """
            models_dir_auto = get_closest_rf_model(dl1, models_dir=models_dir)
            try:
                dl2 = apply_models(dl1, models_dir=models_dir_auto, config=config, telescope=tel, stereo=stereo, mc=ismc)
            except:
                logging.error('RF application failed.')
                os.remove(output_file)
                exit()

            try:
                write_dl2(dl2, output_file=output_file, telescope=tel, config=config, mode='a')
            except:
                logging.error('Writting DL2 file failed.')
                os.remove(output_file)
                exit()

    try:
        write_dl2_info(output_file, rfs_used=models_dir_auto)
    except:
        logging.error('Writting DL2 info failed.')
        os.remove(output_file)
        exit()

if __name__ == '__main__':
    main()