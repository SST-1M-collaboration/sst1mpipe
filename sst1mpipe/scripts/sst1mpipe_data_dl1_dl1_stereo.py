#!/usr/bin/env python

"""
A script to create DL1 stereo data, where each event contains data from both
telescopes. For each DL1 file for tel1, coincident events are searched for
in tel2 DL1 data. Only coincident events are stored in resulting DL1 files.
- Input is a single DL1 file from tel1 and a directory with all relevant 
DL1 files for tel2.
- Output is a DL1 file with coincident events only, to be used in 
sst1mpipe_dl1_dl2 for stereo reconstruction.

Usage:

$> python sst1mpipe_data_dl1_dl1_stereo.py
--input-file-tel1 /tel1/DL1/file_dl1.h5
--input-dir-tel2 /tel2/DL1/
--config sst1mpipe_config.json
--output-dir ./
--pattern dl1_recleaned

"""

import argparse
import glob
import logging
import os
import sys

import sst1mpipe
from sst1mpipe.io import (
    add_wr_dl1_stereo,
    check_outdir,
    load_config,
    load_dl1_sst1m,
    load_more_dl1_tables_mono,
    write_extra_parameters,
)
from sst1mpipe.reco import find_coincidence_offset, make_dl1_stereo
from sst1mpipe.utils import get_pointing_radec, get_stereo_method


def parse_args():
    parser = argparse.ArgumentParser(description="Data DL1 to DL1 stereo")

    # Required arguments
    parser.add_argument(
                        '--input-file-tel1', type=str,
                        dest='input_file_tel1',
                        help='Tel1 DL1 file for which we want to find coincident events in tel2 DL1 data.',
                        required=True
                        )
    parser.add_argument(
                        '--input-dir-tel2', type=str,
                        dest='input_dir_tel2',
                        help='Directory with tel2 DL1 files',
                        required=True
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
                        help='Path to store the output DL1 file with stereo events',
                        default='./'
                        )

    parser.add_argument(
                        '--pattern', type=str,
                        dest='pattern',
                        help='Pattern to select specific files (.e.g. \'dl1_recleaned\') Files from a single night are selected authomaticaly based on the tel 1 filename.',
                        default='dl1'
                        )
    args = parser.parse_args()
    return args

def main():

    args = parse_args()

    input_file_tel1 = args.input_file_tel1
    input_dir_tel2 = args.input_dir_tel2
    outdir = args.outdir
    pattern = args.pattern

    check_outdir(outdir)

    output_file = os.path.join(outdir, input_file_tel1.split('/')[-1].rstrip(".h5") + "_stereo.h5")
    output_logfile = os.path.join(outdir, input_file_tel1.split('/')[-1].rstrip(".h5") + "_stereo.log")

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

    # Expected format of the input DL1 filename is "SST1M{1,2}_{DATE}_*.h5"
    date = input_file_tel1.split('/')[-1].split('_')[1]

    files_selected_t2_0 = sorted(glob.glob(input_dir_tel2 + '/*' + date + '_*_*'+ pattern +'*.h5'))
    files_selected_t2 = [item for item in files_selected_t2_0 if 'pedestal_hist' not in item]

    if len(files_selected_t2) == 0:
        logging.error('No TEL2 file following pattern %s found in %s. Exiting..', pattern, input_dir_tel2)
        exit()
    
    stereo_method = get_stereo_method(config)
    
    if stereo_method == 'SlidingWindow':

        time_offset = find_coincidence_offset(tel1_file=input_file_tel1, 
                                                tel2_files=files_selected_t2, 
                                                outdir=outdir, 
                                                config=config, 
                                                save_figures=True
                                                )
    else: 
        time_offset = 0

    dl1_data_t1 = load_dl1_sst1m(input_file_tel1, tel='tel_021')
    pointing_tel1 = get_pointing_radec(input_file_tel1)

    logging.info('Looking for overlapping tel2 files to the current tel1 file...')
    dl1_data_t2 = load_more_dl1_tables_mono(
        files_selected_t2, 
        config=None, 
        check_finite=False,
        quality_cuts=False,
        time_min=min(dl1_data_t1['local_time']), 
        time_max=max(dl1_data_t1['local_time']),
        pointing_sel=pointing_tel1
        )

    make_dl1_stereo(dl1_file_tel1=input_file_tel1,
                        input_dir_tel2=input_dir_tel2,
                        dl1_data_tel2=dl1_data_t2,
                        file_pattern=pattern,
                        output_path=output_file,
                        config=config,
                        time_offset=time_offset
                        )

    # Write additional params in the DL1 file
    write_extra_parameters(output_file, config=config, ismc=False)

    # Write WR timestamps and pointing back in the DL1 file
    if stereo_method == "WhiteRabbitClosest":
        add_wr_dl1_stereo(output_file, dl1_data_tabs=[dl1_data_t1, dl1_data_t2])

if __name__ == '__main__':
    main()
