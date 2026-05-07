#!/usr/bin/env python

"""
A script to calculate DL1 Hillas Intensity rates per obs ID (~wobble). This is needed 
for run selection, and cannot be done per-file due to statistics.
- Inputs are a directory with DL1 real data files and a directory with DL3 index files.
- Outputs are per-obsid intensity distributions stored in HDF tables.

Usage:

$> python sst1mpipe_extract_dl1_distributions.py
--dl1-dir /data/dl1/
--date 20240120
--output-dir ./
--dl3-index-dir /data/dl3/

"""

import argparse
import glob

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io.misc.hdf5 import write_table_hdf5
from astropy.table import QTable
from gammapy.data import DataStore

from sst1mpipe.io import (
    check_outdir,
    load_dl1_pedestals,
    load_dl1_sst1m,
    load_dl2_sst1m,
)
from sst1mpipe.utils import get_telescopes
from sst1mpipe.utils.NSB_tools import VAR_to_NSB


def parse_args():

    parser = argparse.ArgumentParser(description="Export distributions of DL1 intensity rates per obsid (wobble)")

    parser.add_argument(
                    '--dl1-dir', type=str,
                    dest='dl1_dir',
                    help='Directory with DL1 files.',
                    default='./'
                    )

    parser.add_argument(
                    '--date', type=str,
                    dest='date',
                    help='Date to be processed.',
                    default=None
                    )

    parser.add_argument(
                    '--output-dir', type=str,
                    dest='outdir',
                    help='Path to store the output distributions',
                    default='./'
                    )

    parser.add_argument(
                    '--dl3-index-dir', type=str,
                    dest='dl3_indexes',
                    help='Directory with DL3 indexes.',
                    default=''
                    )

    parser.add_argument(
                '--histogram-bins', type=int,
                dest='bins',
                help='Binninig of the distributions.',
                default=500
                )

    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    outpath = args.outdir
    datastore = args.dl3_indexes
    data_path_dl1 = args.dl1_dir
    nbins = args.bins
    date = args.date

    bins=np.logspace(1, 5, nbins)
    bins_width = bins[1:] - bins[:-1]

    # get dl3
    data_store = DataStore.from_dir(datastore)
    data_store.obs_table.sort('TSTART')

    # get all obsids (wobbles) for given date
    obsids = np.array(data_store.obs_table['OBS_ID']).astype(str)
    obsids_date = [i for i in obsids if date in i]

    # load all DL1 for given date
    files = glob.glob(data_path_dl1 + '/*_stereo.h5')
    stereo = False
    if len(files) > 0:
        stereo = True

    logs = glob.glob(data_path_dl1 + "/*.log")

    if len(files) == 0:
        files = glob.glob(data_path_dl1 + '/*dl1.h5')
    if len(files) == 0:
        files = glob.glob(data_path_dl1 + '/*dl1_recleaned.h5')

    if len(files) == 0:
        files = glob.glob(data_path_dl1 + '/*dl2.h5')

    files.sort()

    # Sometimes first few files are broken
    for file in files:
        tels = get_telescopes(file)
        if len(tels)> 0:
            break

    for tel in tels:
        dl1_data, ped_fractions, recleaned_fractions, times, nsb = load_data(files, logs, config=None, tel=tel, data_level='dl1')

        # NSB calculated for tel2 in stereo is not correct as only tel1 pedetals are propagated to the stereo DL1 files
        if stereo & (str(tel[-1]) == '2'):
            print('NSB level for tel2 cannot be calculed from the stereo DL1 file. The output table is filled with NaNs')
            nsb = np.empty(len(times))
            nsb.fill(np.nan)

        # iterate over obsids for given date and split data
        for obsid in obsids_date:

            row_mask = obsids == obsid
            tmin = data_store.obs_table[row_mask]['TSTART']
            tmax = data_store.obs_table[row_mask]['TSTOP']

            time_mask = (dl1_data['local_time'] > tmin.value[0]) & (dl1_data['local_time'] <= tmax.value[0])
            time_mask_ped = (times > tmin.value[0]) & (times <= tmax.value[0])

            # get intensity distribution per obsid (wobble)

            data, bin_edges = np.histogram(dl1_data[time_mask].camera_frame_hillas_intensity, bins=bins)
            zenith = 90.-np.mean(dl1_data[time_mask].true_alt_tel)

            # elapsed time: sum of the time differences, excluding large ones
            time_diff = dl1_data[time_mask].local_time.diff()
            t_elapsed = np.sum(time_diff[time_diff < 10.]) * u.s
            timestamp = np.mean(dl1_data[time_mask].local_time)

            ped_fraction_mean = np.mean(ped_fractions[time_mask_ped])
            recleaned_fraction_mean = np.mean(recleaned_fractions[time_mask_ped])
            nsb_per_dl3 = np.mean(nsb[time_mask_ped])

            rate = data / t_elapsed
            differential_rate = rate / bins_width

            print(f'Integration time for obsid {obsid}: {t_elapsed.to(u.hour)}')

            t_table = QTable()
            t_table['t_elapsed'] = np.array([t_elapsed.to(u.hour).value])
            z_table = QTable()
            z_table['zenith'] = np.array([zenith])
            tt_table = QTable()
            tt_table['timestamp'] = np.array([timestamp])
            ped_table = QTable()
            ped_table['survived_pedestal_frac'] = np.array([ped_fraction_mean])
            rec_table = QTable()
            rec_table['recleaned_fraction'] = np.array([recleaned_fraction_mean])
            nsb_table = QTable()
            nsb_table['NSB'] = np.array([nsb_per_dl3])

            # Save all to output tables
            table = get_table(rate, differential_rate, bin_edges=bin_edges)

            check_outdir(outpath)
            outfile = outpath + '/intensity_hist_'+tel+'_'+str(int(obsid))+'.h5'
            write_table_hdf5(table, outfile, path='intensity_hist', overwrite=True, append=True, serialize_meta=True)
            write_table_hdf5(t_table, outfile, path='t_elapsed', overwrite=True, append=True, serialize_meta=True)
            write_table_hdf5(z_table, outfile, path='zenith', overwrite=True, append=True, serialize_meta=True)
            write_table_hdf5(tt_table, outfile, path='timestamp', overwrite=True, append=True, serialize_meta=True)
            write_table_hdf5(ped_table, outfile, path='survived_pedestal_frac', overwrite=True, append=True, serialize_meta=True)
            write_table_hdf5(rec_table, outfile, path='recleaned_fraction', overwrite=True, append=True, serialize_meta=True)
            write_table_hdf5(nsb_table, outfile, path='NSB', overwrite=True, append=True, serialize_meta=True)

def get_table(rate, diff_rate, bin_edges=None):

    table = QTable()
    table["low"] = bin_edges[:-1]
    table["high"] = bin_edges[1:]
    table["center"] = (bin_edges[:-1] + bin_edges[1:]) / 2
    table["rate"] = rate
    table["diff_rate"] = diff_rate

    return table

def load_data(files, logs, config=None, tel=None, data_level='dl1'):
    i = 0
    data = None
    ped_fractions = []
    recleaned_fractions = []
    times = []
    nsb = []

    for input_file in files:
        try:
            if data_level == 'dl1':
                df = load_dl1_sst1m(input_file, tel=tel, config=config, table='pandas')
                times.append(df['local_time'][0])

                # find log
                date = input_file.split('/')[-1].split('.')[0].split('_')[1]
                run = input_file.split('/')[-1].split('.')[0].split('_')[2]
                res = [i for i in logs if date in i]
                log_file = [i for i in res if run in i][0]

                # find fraction of pedestals in log
                ped = 'Fraction of pedestal'
                reclean = 'raised picture threshold'
                with open(log_file) as fp:
                    # read all lines in a list
                    lines = fp.readlines()
                    ped_fraction = 100
                    recleaned_fraction = 0
                    for line in lines:
                        # check if string present on a current line
                        if line.find(ped) != -1:
                            ped_fraction = float(line.split(': ')[-1])
                        if line.find(reclean) != -1:
                            recleaned_fraction = float(line.split(': ')[-1])
                    ped_fractions.append(ped_fraction)
                    recleaned_fractions.append(recleaned_fraction)

            elif data_level == 'dl2':
                df = load_dl2_sst1m(input_file, tel=tel, config=config, table='pandas')
        except Exception:
            print('Broken file: ' + input_file + ', skipping.')
            continue
        try:
            pt = load_dl1_pedestals(input_file)
            if '1' in tel:
                cs=21
            else:
                cs=22
            NSB = VAR_to_NSB(pt['pedestal_charge_std'].mean(axis=1)**2, cs)
            nsb.append(NSB.mean())
        except Exception:
            print('No pedestals in : ' + input_file + '.')
            nsb.append(np.nan)

        if i == 0:
            data = df
            i += 1
        else:
            data = pd.concat([data, df])

    ped_fractions = np.array(ped_fractions)
    recleaned_fractions = np.array(recleaned_fractions)
    times = np.array(times)
    nsb = np.array(nsb)
    return data, ped_fractions, recleaned_fractions, times, nsb

if __name__ == '__main__':
    main()


