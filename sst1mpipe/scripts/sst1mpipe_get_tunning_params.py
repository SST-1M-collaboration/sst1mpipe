#!/usr/bin/env python

"""
A script to calculate NSB tunning parameters for MC.
- Inputs are charge distribution files calculated from pedestal events in real data 
(usualy in the same dir as DL1 tables), and a directory with charge distributions
from MC.
- Outputs are tunning parameters for each data file, and images comparing MC/data 
charge distributions.

Usage:

$> python sst1mpipe_get_tunning_params.py
--input-files /data/dl1/*hist.h5
--mc-hist-dir /mc/dl1/px_charges/
--pattern mc_hist*.h5
--output-dir ./

"""

import argparse
import glob

import matplotlib.pyplot as plt
import numpy as np

from sst1mpipe.io import check_outdir, read_charges_data, read_charges_mc


def parse_args():
    parser = argparse.ArgumentParser(description="Extraction of MC tunning parameters from pixel charge histograms")

    # Required arguments
    parser.add_argument(
                        '--input-files', '-f', type=str,
                        dest='input_files',
                        help='Input histogram files from DATA',
                        required=True,
                        nargs='*'
                        )

    parser.add_argument(
                        '--mc-hist-dir', '-m', type=str,
                        dest='mc_hist_dir',
                        help='Directory with MC histograms (to be loaded all due to low statistics in individual files)',
                        required=True
                        )

    parser.add_argument(
                        '--pattern', type=str,
                        dest='pattern',
                        help='File pattern for MC pixel histograms',
                        default='*_hist.h5'
                        )

    parser.add_argument(
                        '--output-dir', '-o', type=str,
                        dest='outdir',
                        help='Path to store the output DL1 file',
                        default='./'
                        )

    args = parser.parse_args()
    return args

def main():

    args = parse_args()
    input_data_hists = args.input_files
    mc_hist_dir = args.mc_hist_dir
    pattern = args.pattern
    outdir = args.outdir

    check_outdir(outdir)

    mc_hist_files = glob.glob(mc_hist_dir + '/' + pattern)
    
    i = 0
    for mc_hist_file in mc_hist_files:

        bin_centers_mc, charge_mc = read_charges_mc(mc_hist_file, tel=1)
        if i == 0:
            i += 1
            charges_mc = charge_mc
        else:
            charges_mc += charge_mc

    charges_mc_norm = charges_mc / sum(charges_mc)
    mode_mc = bin_centers_mc[np.argmax(charges_mc_norm)]

    biases, noises, file_numbers = [], [], []

    for input_data_hist in input_data_hists:

        print(input_data_hist)
        bin_centers, charge = read_charges_data(input_data_hist)

        # normalization of histograms
        charges_data_norm = charge / sum(charge)
        
        # determination of the shift of the distribution (bias)
        mode_data = bin_centers[np.argmax(charges_data_norm)]

        bias = mode_data - mode_mc
        bias = max(0., bias)
        biases.append(bias)

        # determination of widening of the nsb bump
        data_diff = charges_data_norm - mode_data
        mc_diff = charges_mc_norm - mode_mc

        # cutting out outliers in data and shower pixels in MC
        max_diff = 6
        # calculate widening of the noise bump:
        added_noise = (np.sum(data_diff[abs(data_diff)<max_diff]**2)/len(data_diff[abs(data_diff)<max_diff]) 
                        - np.sum(mc_diff[abs(mc_diff)<max_diff]**2)/len(mc_diff[abs(mc_diff) < max_diff]))
        noise = max(0., added_noise)
        noises.append(noise)

        file_no = int(input_data_hist.split('/')[-1].split('_')[2])
        file_numbers.append(file_no)

        print(f'File: {file_no}, Bias: {bias}, Added Noise: {noise}')

        plt.figure()
        plt.plot(bin_centers, charges_data_norm, label='data')
        plt.plot(bin_centers_mc, charges_mc_norm, label='mc')
        plt.axvline(mode_data)
        plt.yscale('log')
        plt.xlabel('pixel charge [p.e.]')
        plt.xlim([-2, 10])
        plt.grid()
        plt.legend()
        outfile = input_data_hist.split('/')[-1].split('.')[0] + '.png'
        plt.savefig(outdir + '/' + outfile)
        plt.close()

    plt.figure()
    plt.plot(file_numbers, noises, label='added noise')
    plt.plot(file_numbers, biases, label='bias')
    plt.xlabel('File number')
    plt.ylabel('pixel charge [p.e.]')
    plt.grid()
    plt.legend()
    plt.xlabel('File')
    plt.savefig(outdir + '/all_params.png')
    plt.close()

if __name__ == '__main__':
    main()
