#!/usr/bin/env python


import sst1mpipe
import argparse
import sys
import os
import logging
from sst1mpipe.io import (
    load_dl2_sst1m,
    load_config,
    check_outdir,
    write_dl2
)

from sst1mpipe.utils import (
    get_telescopes, 
    check_mc
)
from sst1mpipe.reco import (
    reco_misdirection,
    get_evttype_edges,
    classify_evt_types
)

import numpy as np
from astropy.io.misc.hdf5 import write_table_hdf5, read_table_hdf5
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description="MC train RFs")

    # Required arguments
    parser.add_argument(
                        '--input-file', '-f',
                        dest='file',
                        required=True,
                        help='Path to any DL2 file. On data, the script reconstructs misdirection and classify events. On MC, it can do the same (useful for IRFs), or it calculates edges for event classification.'
                        )
    parser.add_argument('--config', '-c', action='store', type=str,
                        dest='config_file',
                        help='Path to a configuration file.',
                        required=True
                        )
    parser.add_argument(
                        '--mis-models-dir',
                        action='store',
                        help='Path to models for reconstruction of misdirection.',
                        dest='misdirection_models',
                        required=True
                        )
    # Optional arguments
    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store the edges for classification.',
                        default='./'
                        )
    parser.add_argument(
                        '--mis-edges-dir',
                        action='store',
                        help='Directory with stored misdirection edges to classify the events. If NOT used, the script will calculate energy dependent edges from distribution of reco misdirection and store them in --output-dir.',
                        dest='edges_dir',
                        default=None
                        )

    args = parser.parse_args()
    return args

def main():

    args = parse_args()

    input_file = args.file
    outdir = args.outdir
    misdirection_models_dir = args.misdirection_models
    edges_dir = args.edges_dir

    check_outdir(outdir)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers= [
            logging.FileHandler(outdir+'/sst1mpipe_event_classes.log', 'w+'),
            logging.StreamHandler(stream=sys.stdout)
            ]
    )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    config = load_config(args.config_file)

    tels = get_telescopes(input_file, data_level="dl2")
    ismc = check_mc(input_file)

    for tel in tels:
        
        dl2 = load_dl2_sst1m(input_file, tel=tel, config=config, table='pandas')
        dl2 = reco_misdirection(dl2, models_dir=misdirection_models_dir, config=config, telescope=tel)

        if ismc and edges_dir is None:
            logging.info('Determining edges of misdirection classes..')
            edges = get_evttype_edges(dl2, config=config, percentiles=[25, 50, 75])
            outfile = outdir + '/misdirection_edges_'+tel+'.h5'
            write_table_hdf5(edges, outfile, path='misdirection_edges', overwrite=True, append=True, serialize_meta=True)
            logging.info('Done. Edges stored in '+outfile)

            dl2 = classify_evt_types(dl2, edges=edges)

            evttypes = dl2['event_type'].unique()
            fig, ax = plt.subplots(1, 1, figsize=(7, 6))
            for etype in evttypes:
                mask = dl2['event_type'] == etype
                ax.plot(dl2['reco_energy'][mask], dl2['log_reco_misdirection'][mask], '.', alpha=0.2, label='Evt type ' + str(etype))
            ax.set_xscale('log')
            ax.set_xlabel('Reco energy [TeV]')
            ax.set_ylabel('log Reco misdirection [deg]')
            ax.legend()
            plt.tight_layout()
            plt.savefig(outdir + 'classes_'+tel+'.png', dpi=200)

        elif edges_dir is not None:
            logging.info('Determining event classes in the file.')

            edges_file = edges_dir + '/misdirection_edges_'+tel+'.h5'
            edges = read_table_hdf5(edges_file, path='misdirection_edges')
            dl2 = classify_evt_types(dl2, edges=edges)
            # We just update the DL2 table in the original file with new columns
            write_dl2(dl2, output_file=input_file, telescope=tel, config=config, mode='a')

        else:
            logging.info('Invalid combination of MC/data input file and --get-edges parameter. See --help.')

        
if __name__ == '__main__':
    main()