#!/usr/bin/env python

"""
A script to create night summary for each observed source, which is 
stored as a PDF in the base data analysis directory. For the full
summary, the script needs all data analysis levels up to DL3 for 
mono and stereo to be already produced. It also needs the DL1 
distributions to be already extracted.
- Input is a base data analysis directory, within which the expected 
structure is BASE_DIR/{date}/{source}/{telescope}/{data level}/{subdir}/
- Outputs are night summary PDFs for all sources observed during the night

Usage:

$> python sst1mpipe_night_summary.py
--base-dir /data/
--out-dir /outputs/
--date 20240810
--config sst1mpipe_config.json
--sub-dir /v0.5.5/
--source-catalog sst1m_source_catalog.json

"""

import argparse
import sys
import os
import logging

import sst1mpipe

from sst1mpipe.io import (
    load_dl1_sst1m,
    load_dl2_sst1m,
    load_dl1_pedestals,
    load_config,
    load_source_catalog,
    check_outdir
)
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sst1mpipe.utils import (
    get_moon_params,
    get_wr_timestamp,
    get_sources_in_dir
)

from sst1mpipe.io import load_distributions_sst1m

import astropy.units as u
from astropy.time import Time
from sst1mpipe.utils.NSB_tools import plot_average_nsb_VS_time
from astropy.table import vstack
from PIL import Image
from gammapy.maps import MapAxis

from sst1mpipe.analysis import get_theta2_from_dl3, plot_theta2_dl3
from astropy.coordinates import SkyCoord


def parse_args():
    parser = argparse.ArgumentParser(description="Makes run summary PDF for all sources observed during a single night.")

    parser.add_argument(
                        '--base-dir', '-d',
                        dest='base_dir',
                        required=True,
                        type=str,
                        help='Base data analysis directory.'
                        )
    parser.add_argument(
                        '--out-dir', '-o',
                        dest='out_dir',
                        type=str,
                        help='Output directory to store the night summary.',
                        default=''
                        )

    parser.add_argument(
                        '--date',
                        dest='date',
                        default=None,
                        type=str,
                        help='Date of the observing night (YYYYMMDD).',
                        required=True
                        )

    parser.add_argument('--config', '-c', action='store', type=str,
                        dest='config_file',
                        help='Path to a configuration file.',
                        required=True
                        )

    parser.add_argument(
                        '--sub-dir', '-s', type=str,
                        dest='version',
                        help='Sub directory for given data sample inside each data level directory, e.g. a version of sst1mpipe used (../DL{1,2,3}/v0.5.5/).',
                        default=''
                        )

    parser.add_argument('--source-catalog', action='store', type=str,
                        dest='source_cat',
                        help='Path to a json file with catalog of all observed SST1M targets.',
                        default=''
                        )

    parser.add_argument(
                        '--read-dl1-from-dl2',
                        action='store_true',
                        help='If used, the script will ignore the DL1 directories and tries to read DL1 tabs from DL2 files.',
                        dest='read_dl1_tab_from_dl2_files'
                        )

    args = parser.parse_args()
    return args


def load_files(files, config=None, tel=None, level='dl1', stereo=False):
    i = 0
    ped_table = []
    data = []

    for input_file in files:
        try:
            if level == 'dl1':
                df = load_dl1_sst1m(input_file, tel=tel, config=config, table='pandas', stereo=stereo)
                if not stereo:
                    try:
                        pt = load_dl1_pedestals(input_file)
                    except: 
                        pass
            elif level == 'dl2':
                df = load_dl2_sst1m(input_file, tel=tel, config=config, table='pandas')
        except:
            logging.warning("Broken file %s", input_file)
            continue
        if i == 0:
            data = df
            if (level == 'dl1') and not stereo:
                try:
                    ped_table = pt
                except:
                    pass
        else:
            try:
                data = pd.concat([data, df])
            except:
                logging.warning("Broken file %s", input_file)
                continue
            if (level == 'dl1') and not stereo:
                try:
                    ped_table = vstack([ped_table, pt])
                except:
                    logging.warning("No pedestal monitoring in file %s", input_file)
                    continue
        i += 1
    return data, ped_table


def make_file_list(path, stereo=False, level='dl1'):
    
    files_selected = []
    
    if stereo:
        if level == 'dl1':
            files = glob.glob(path + '/*dl1_stereo.h5')
        else:
            files = glob.glob(path + '/*stereo_dl2.h5')
    else:
        files = glob.glob(path + '/*'+level+'.h5')

    files_selected = files_selected + files
    
    files_selected = np.sort(files_selected)
    return files_selected


def prepare_to_pdf(image_path):

    img1 = Image.open(image_path)
    img1.load() # required for png.split()
    background1 = Image.new("RGB", img1.size, (255, 255, 255))
    background1.paste(img1, mask=img1.split()[3])

    return background1


def get_min_max_times(dl1_files, tel=None, config=None, stereo=False):

    dl1_files = np.sort(dl1_files)
    ind = 0
    first_loaded = False
    while not first_loaded:
        try:
            df = load_dl1_sst1m(dl1_files[ind], tel=tel, config=config, table='pandas', stereo=stereo)
            first_loaded = True
        except:
            logging.warning('Broken file')
            ind += 1
            if ind > 50: 
                first_loaded = True
                logging.warning('Something is really wrong, too many broken files!')
    min_time = min(df['local_time'])
    ind = 0
    last_loaded = False
    while not last_loaded:
        try:
            df = load_dl1_sst1m(dl1_files[len(dl1_files)-ind], tel=tel, config=config, table='pandas', stereo=stereo)
            last_loaded = True
        except:
            logging.warning('Broken file')
            ind += 1
            if ind > 50: 
                last_loaded = True
                logging.warning('Something is really wrong, too many broken files!')
    max_time = max(df['local_time'])
    return min_time, max_time


def main():

    args = parse_args()

    base_path = args.base_dir
    date = args.date
    config_file = args.config_file
    version = args.version
    out_dir = args.out_dir
    source_catalog_file = args.source_cat
    read_dl1_tab_from_dl2_files = args.read_dl1_tab_from_dl2_files

    if read_dl1_tab_from_dl2_files:
        dl1_tab_data_level = 'dl2'
        logging.warning('DL1 tables read from DL2 files.')
    else:
        dl1_tab_data_level = 'dl1'


    if len(out_dir) > 0:
        outpath = out_dir
    else:
        outpath = base_path + '/' + date + '/'
        print('Default outpath used: %s', outpath)

    check_outdir(outpath)

    logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers= [
                logging.FileHandler(outpath+'/sst1mpipe_night_summary_'+str(sst1mpipe.__version__) + '.log', 'w+'),
                logging.StreamHandler(stream=sys.stdout)
                ]
        )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    sources = get_sources_in_dir(base_path + '/' + date + '/')

    config = load_config(config_file)

    is_catalog = False
    if len(source_catalog_file):
        try:
            source_catalog = load_source_catalog(source_catalog_file)
            is_catalog = True
        except:
            logging.warning('Source catalog file not found!')
    else:
        logging.warning('Source catalog file not specified!')

    tels = ['cs1', 'cs2', 'stereo']

    for source in sources:

        logging.info(' ==== Preparing summary for source %s ==== ', source)
        fig, ax = plt.subplots(1, 1, figsize=(12, 5)) # dl1 event rates
        fig1, ax1 = plt.subplots(1, 1, figsize=(12, 5)) # dl1 event rates (zoom)
        fig2, ax2 = plt.subplots(1, 1, figsize=(12, 5)) # dl2 event rates
        fig3, ax3 = plt.subplots(1, 1, figsize=(12, 5)) # NSB
        fig4, ax4 = plt.subplots(1, 2, figsize=(12,5)) # distributions, tel1
        fig5, ax5 = plt.subplots(1, 2, figsize=(12,5)) # distributions, tel2
        fig6, ax6 = plt.subplots(1, 2, figsize=(12,5)) # distributions, stereo
        fig7, ax7 = plt.subplots(1, 2, figsize=(12, 5)) # zenith angles
        fig8, ax8 = plt.subplots(1, 2, figsize=(12, 5)) # cogs, dl1 mono, survived cleaning
        fig9, ax9 = plt.subplots(1, 2, figsize=(12, 5)) # cogs, dl2 mono, photons
        fig10, ax10 = plt.subplots(1, 2, figsize=(12, 5)) # DL1 stereo, survived cleaning
        fig11, ax11 = plt.subplots(1, 1, figsize=(12, 5)) # Moon
        fig12, ax12 = plt.subplots(1, 2, figsize=(12, 5)) # stereo dt
        fig13, ax13 = plt.subplots(1, 3, figsize=(16, 5)) # theta2 plots

        median1 = []
        median2 = []
        t_diff_all = []
        times_all = []
        moon_done = False

        rate_colors = {'tel_021':'blue', 'tel_022':'orange', 'stereo':'green'}

        logging.info('Making DL1 time bins to speed up the processing.')
        bunch_size = 20 # N of files

        for telescope in ['cs1', 'cs2']:
            logging.info('Telescope %s', telescope)
            base_path_source = base_path + '/' + date + '/' + source+'/' + telescope
            dl1_path = base_path_source + '/' + dl1_tab_data_level.upper() + '/' + version
            is_dl1 = len(glob.glob(dl1_path+'/'+'*.h5'))
            # handle lower case directories
            if not is_dl1:
                dl1_path = base_path_source + '/' + dl1_tab_data_level.lower() + '/' + version
                is_dl1 = len(glob.glob(dl1_path+'/'+'*.h5'))
                
            if telescope == 'cs1':
                if is_dl1:
                    dl1_files = make_file_list(dl1_path, stereo=False, level=dl1_tab_data_level)
                    tel1_min, tel1_max = get_min_max_times(dl1_files, tel='tel_021', config=config)
                else:
                    tel1_min = 1e15
                    tel1_max = 0
            elif telescope == 'cs2':
                if is_dl1:
                    dl1_files = make_file_list(dl1_path, stereo=False, level=dl1_tab_data_level)
                    tel2_min, tel2_max = get_min_max_times(dl1_files, tel='tel_022', config=config)
                else:
                    tel2_min = 1e15
                    tel2_max = 0
        min_time_dl1 = min([tel1_min, tel2_min])
        max_time_dl1 = max([tel1_max, tel2_max])
        if max_time_dl1 == 0:
            logging.error('Determination of time bins from %s files failed.', dl1_tab_data_level.upper())
            exit()
        dl1_rate_bins = np.linspace(min_time_dl1, max_time_dl1, int((max_time_dl1-min_time_dl1)/10.))
        logging.info('Min_time %f, Max time %f', min_time_dl1, max_time_dl1)

        for telescope in tels:

            # Load all files
            logging.info(' ==== Preparing summary for telescope %s ==== ', telescope)
            base_path_source = base_path + '/' + date + '/' + source+'/' + telescope
            dl1_path = base_path_source + '/'+dl1_tab_data_level.upper()+'/' + version
            dl2_path = base_path_source + '/DL2/' + version
            dl3_path = base_path_source + '/DL3/' + version
            dist_path = base_path_source + '/distributions/' + version
            
            is_dl1 = len(glob.glob(dl1_path+'/'+'*.h5'))
            is_dl2 = len(glob.glob(dl2_path+'/'+'*.h5'))
            is_dl3 = len(glob.glob(dl3_path+'/'+'*.fits'))
            is_dist = len(glob.glob(dist_path+'/'+'*.h5'))

            # handle low case directories
            if not is_dl1:
                dl1_path = base_path_source + '/'+dl1_tab_data_level.lower()+'/' + version
                is_dl1 = len(glob.glob(dl1_path+'/'+'*.h5'))
            if not is_dl2:
                dl2_path = base_path_source + '/dl2/' + version
                is_dl2 = len(glob.glob(dl2_path+'/'+'*.h5'))
            if not is_dl3:
                dl3_path = base_path_source + '/dl3/' + version
                is_dl3 = len(glob.glob(dl3_path+'/'+'*.fits'))

            if not is_dl1:
                logging.warning('DL1 data is missing, there will be some missing figures in the final PDF.')
            if not is_dl2:
                logging.warning('DL2 data is missing, there will be some missing figures in the final PDF.')
            if not is_dl3:
                logging.warning('DL3 data is missing, there will be some missing figures in the final PDF.')
            if not is_dist:
                logging.warning('DL1 rate distributions are missing, there will be some missing figures in the final PDF.')
            
            if telescope not in 'stereo':
                stereo=False
            else:
                stereo=True
            
            if is_dl1:
                dl1_files = make_file_list(dl1_path, stereo=stereo, level=dl1_tab_data_level)
            if is_dl2:
                dl2_files = make_file_list(dl2_path, stereo=stereo, level='dl2')
            
            if '1' in telescope:
                tel = 'tel_021'
                tt = 21
            elif '2' in telescope:
                tel = 'tel_022'
                tt = 22
            else:
                tel = 'stereo'
                tt = 0

            if is_dl1:
                # We load and plot DL1 rates in bunches of file, 
                # because otherwise it needs too much memory and
                # is realy slow.
                N_bunches = int(len(dl1_files)/bunch_size)
                h11_tot = np.zeros((100, 100))
                h22_tot = np.zeros((100, 100))
                if tel != 'stereo':
                    h_tot = np.zeros(len(dl1_rate_bins)-1)
                    for i in range(N_bunches):
                        logging.info('BUNCH %d', i)
                        dl1, ped_table = load_files(dl1_files[i*bunch_size:bunch_size*(i+1)], tel=tel, level='dl1')
                        # Trigger rates
                        if len(dl1) > 0:
                            h1 = np.histogram(dl1.local_time, bins=dl1_rate_bins)
                            h_tot += h1[0]
                            # pedestal table
                            if len(ped_table) > 0:
                                plot_average_nsb_VS_time(ped_table,tt,ax=ax3, color=rate_colors[tel], label=None)
                            # CoGs - DL1 mono, survived cleaning
                            if tel == 'tel_021':
                                h11, xedges, yedges = np.histogram2d(dl1['camera_frame_hillas_x'].dropna(), dl1['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                                h11_tot += h11
                            else:
                                h22, xedges, yedges = np.histogram2d(dl1['camera_frame_hillas_x'].dropna(), dl1['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                                h22_tot += h22
                    # last bunch
                    dl1, ped_table = load_files(dl1_files[bunch_size*(i+1):], tel=tel, level='dl1')
                    if len(dl1) > 0:
                        h1 = np.histogram(dl1.local_time, bins=dl1_rate_bins)
                        h_tot += h1[0]
                    centers = (dl1_rate_bins[1:]+dl1_rate_bins[:-1]) / 2
                    ax.plot(centers, h_tot/10, label=tel, alpha=0.7, color=rate_colors[tel])
                    ax1.plot(centers, h_tot/10, label=tel, alpha=0.7, color=rate_colors[tel])
                    median1.append(np.median(h_tot/10.))
                    # pedestal table
                    if len(ped_table) > 0:
                        plot_average_nsb_VS_time(ped_table,tt,ax=ax3, color=rate_colors[tel], label=tel)
                        # CoGs - DL1 mono, survived cleaning
                    X, Y = np.meshgrid(xedges, yedges)
                    if tel == 'tel_021':
                        if len(dl1) > 0:
                            h11, xedges, yedges = np.histogram2d(dl1['camera_frame_hillas_x'].dropna(), dl1['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                            h11_tot += h11
                        ax8[0].pcolormesh(X, Y, h11_tot)
                    else:
                        if len(dl1) > 0:
                            h22, xedges, yedges = np.histogram2d(dl1['camera_frame_hillas_x'].dropna(), dl1['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                            h22_tot += h22
                        ax8[1].pcolormesh(X, Y, h22_tot)
                else:
                    h_tot = np.zeros(len(dl1_rate_bins)-1)
                    for i in range(N_bunches):
                        dl1_stereo_1,_ = load_files(dl1_files[i*bunch_size:bunch_size*(i+1)], tel='tel_021', level='dl1', stereo=stereo)
                        dl1_stereo_2,_ = load_files(dl1_files[i*bunch_size:bunch_size*(i+1)], tel='tel_022', level='dl1', stereo=stereo)
                        # Trigger rates
                        if (len(dl1_stereo_1) > 0) and (len(dl1_stereo_1) == len(dl1_stereo_2)):
                            h1 = np.histogram(dl1_stereo_1.local_time, bins=dl1_rate_bins)
                            h_tot += h1[0]
                            # CoG - DL1 stereo, survived cleaning
                            h11, xedges, yedges = np.histogram2d(dl1_stereo_1['camera_frame_hillas_x'].dropna(), dl1_stereo_1['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                            h22, xedges, yedges = np.histogram2d(dl1_stereo_2['camera_frame_hillas_x'].dropna(), dl1_stereo_2['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                            h11_tot += h11
                            h22_tot += h22
                            # WR timestamps
                            t_t1 = get_wr_timestamp(dl1_stereo_1)
                            t_t2 = get_wr_timestamp(dl1_stereo_2)
                            dt = t_t1 - t_t2
                            t_diff_all.append(dt)
                            times_all.append(dl1_stereo_1.local_time.to_numpy())
                            #print(t_diff_all)
                            #print(len(times_all))

                    # last bunch
                    dl1_stereo_1,_ = load_files(dl1_files[bunch_size*(i+1):], tel='tel_021', level='dl1', stereo=stereo)
                    dl1_stereo_2,_ = load_files(dl1_files[bunch_size*(i+1):], tel='tel_022', level='dl1', stereo=stereo)
                    if (len(dl1_stereo_1) > 0) and (len(dl1_stereo_1) == len(dl1_stereo_2)):
                        h1 = np.histogram(dl1_stereo_1.local_time, bins=dl1_rate_bins)
                        h_tot += h1[0]
                        h11, xedges, yedges = np.histogram2d(dl1_stereo_1['camera_frame_hillas_x'].dropna(), dl1_stereo_1['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                        h22, xedges, yedges = np.histogram2d(dl1_stereo_2['camera_frame_hillas_x'].dropna(), dl1_stereo_2['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                        h11_tot += h11
                        h22_tot += h22
                        # WR timestamps
                        t_t1 = get_wr_timestamp(dl1_stereo_1)
                        t_t2 = get_wr_timestamp(dl1_stereo_2)
                        dt = t_t1 - t_t2
                        t_diff_all.append(dt)
                        times_all.append(dl1_stereo_1.local_time.to_numpy())

                    centers = (dl1_rate_bins[1:]+dl1_rate_bins[:-1]) / 2
                    ax.plot(centers, h_tot/10, label=tel, alpha=0.7, color=rate_colors[tel])
                    ax1.plot(centers, h_tot/10, label=tel, alpha=0.7, color=rate_colors[tel])
                    median1.append(np.median(h_tot/10.))
                    X, Y = np.meshgrid(xedges, yedges)
                    ax10[0].pcolormesh(X, Y, h11_tot)
                    ax10[1].pcolormesh(X, Y, h22_tot)
                    t_diff_all = np.concatenate(t_diff_all)
                    times_all = np.concatenate(times_all)
                    ax12[0].hist(t_diff_all, bins=100, range=[-1000, 1000])
                    ax12[1].plot(times_all, t_diff_all, '.')

            if is_dl2:
                # Reconstructed event rates (zoomed)
                h_tot = np.zeros(len(dl1_rate_bins)-1)
                h11_tot = np.zeros((100, 100))
                h22_tot = np.zeros((100, 100))
                zeniths_all = np.zeros(50)
                local_time = []
                zenith_time = []
                time_all = []
                moon_altaz_all = []
                moon_separation_all = []
                moon_phase_angle_all = []
                N_bunches = int(len(dl2_files)/bunch_size)
                for i in range(N_bunches):
                    logging.info('BUNCH %d', i)
                    dl2,_ = load_files(dl2_files[i*bunch_size:bunch_size*(i+1)], tel=tel, level='dl2')
                    if len(dl2) > 0:
                        h1 = np.histogram(dl2.local_time, bins=dl1_rate_bins)
                        h_tot += h1[0]

                    # Zenith angles
                    if (tel != 'stereo') and (len(dl2) > 0):
                        zenith = 90-dl2['true_alt_tel']
                        h1 = np.histogram(zenith, bins=50, range=[0, 65])
                        zeniths_all += h1[0]
                        local_time.append(np.array(dl2['local_time'].T))
                        zenith_time.append(np.array(zenith.T))

                    # moon
                    if (not moon_done) and (len(dl2) > 0):
                        time, moon_altaz, moon_separation, moon_phase_angle = get_moon_params(dl2, config=config, tel=tel, thinning=100)
                        time_all.append(np.array(time.unix))
                        moon_altaz_all.append(np.array(moon_altaz.alt.to_value(u.deg)))
                        moon_separation_all.append(np.array(moon_separation.to_value(u.deg)))
                        moon_phase_angle_all.append(np.array(moon_phase_angle.to_value(u.deg)))

                    # CoGs - DL2 mono, gammas
                    if len(dl2) > 0:
                        mask1 = dl2['gammaness'] > 0.7
                        if tel == 'tel_021':
                            h11, xedges, yedges = np.histogram2d(dl2[mask1]['camera_frame_hillas_x'].dropna(), dl2[mask1]['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                            h11_tot += h11
                        elif tel == 'tel_022':
                            h22, xedges, yedges = np.histogram2d(dl2[mask1]['camera_frame_hillas_x'].dropna(), dl2[mask1]['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                            h22_tot += h22

                # last bunch
                dl2,_ = load_files(dl2_files[bunch_size*(i+1):], tel=tel, level='dl2')
                if len(dl2) > 0:
                    h1 = np.histogram(dl2.local_time, bins=dl1_rate_bins)
                    h_tot += h1[0]
                centers = (dl1_rate_bins[1:]+dl1_rate_bins[:-1]) / 2
                ax2.plot(centers, h_tot/10., alpha=0.7, label=tel, color=rate_colors[tel])
                median2.append(np.median(h_tot/10.))
                if tel != 'stereo':
                    if len(dl2) > 0:
                        zenith = 90-dl2['true_alt_tel']
                        h1 = np.histogram(zenith, bins=50, range=[0, 65])
                        zeniths_all += h1[0]
                        local_time.append(np.array(dl2['local_time'].T))
                        zenith_time.append(np.array(zenith.T))
                    ax7[0].step(np.linspace(0, 65, 50), zeniths_all, alpha=0.5, label=tel, where='mid')
                    local_time = np.concatenate(local_time)
                    zenith_time = np.concatenate(zenith_time)
                    ax7[1].plot(local_time, zenith_time, '.', label=tel, alpha=0.7)
                if not moon_done:
                    if len(dl2) > 0:
                        time, moon_altaz, moon_separation, moon_phase_angle = get_moon_params(dl2, config=config, tel=tel, thinning=100)
                        time_all.append(np.array(time.unix))
                        moon_altaz_all.append(np.array(moon_altaz.alt.to_value(u.deg)))
                        moon_separation_all.append(np.array(moon_separation.to_value(u.deg)))
                        moon_phase_angle_all.append(np.array(moon_phase_angle.to_value(u.deg)))

                    ax11.plot(np.concatenate(time_all), np.concatenate(moon_altaz_all), label='Moon alt')
                    ax11.plot(np.concatenate(time_all), np.concatenate(moon_separation_all), label='Moon sep')
                    ax11.plot(np.concatenate(time_all), np.concatenate(moon_phase_angle_all), label='Moon phase (full=0)')
                    moon_done = True

                # CoGs - DL2 mono, gammas
                X, Y = np.meshgrid(xedges, yedges)
                if len(dl2) > 0:
                    mask1 = dl2['gammaness'] > 0.7
                if tel == 'tel_021':
                    if len(dl2) > 0:
                        h11, xedges, yedges = np.histogram2d(dl2[mask1]['camera_frame_hillas_x'].dropna(), dl2[mask1]['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                        h11_tot += h11
                    ax9[0].pcolormesh(X, Y, h11_tot)
                elif tel == 'tel_022':
                    if len(dl2) > 0:
                        h22, xedges, yedges = np.histogram2d(dl2[mask1]['camera_frame_hillas_x'].dropna(), dl2[mask1]['camera_frame_hillas_y'].dropna(), bins=100, range=[[-0.5, 0.5], [-0.5, 0.5]])
                        h22_tot += h22
                    ax9[1].pcolormesh(X, Y, h22_tot)

            # DL3 Datastore
            plot_theta2 = False
            if is_dl3:
                try:
                    target_coords = SkyCoord.from_name(source)
                    plot_theta2 = True
                except:
                    logging.info('%s coordinates cannot be guessed authomaticaly, using source catalog file %s.', source, source_catalog_file)
                    if is_catalog:
                        if len(source_catalog[source]):
                            logging.info('Source coordinates found in the catalog file: %s', source_catalog[source])
                            if source_catalog[source]['frame'] == 'icrs':
                                target_coords = SkyCoord(
                                    ra=source_catalog[source]['ra']*u.deg, 
                                    dec=source_catalog[source]['dec']*u.deg, 
                                    frame='icrs'
                                    )
                                logging.info('ra: %f, dec: %f', source_catalog[source]['ra'], source_catalog[source]['dec'])
                                plot_theta2 = True
                            elif source_catalog[source]['frame'] == 'galactic':
                                target_coords = SkyCoord(
                                    l=source_catalog[source]['l']*u.deg, 
                                    b=source_catalog[source]['b']*u.deg, 
                                    frame='galactic'
                                    )
                                logging.info('l: %f, b: %f', source_catalog[source]['l'], source_catalog[source]['b'])
                                plot_theta2 = True
                        else:
                            logging.warning('Source %s not found in the source catalog, unknown source coords, skipping theta2 figures.', source)
                    else:
                        logging.warning('Source catalog file %s not found, unknown source coords, skipping theta2 figures.', source_catalog_file)
                if plot_theta2:
                    theta2_axis = MapAxis.from_bounds(0, 0.5, nbin=10, interp="lin", unit="deg2")
                    theta_cut = 0.2 * u.deg
                    counts_on, counts_off, alpha, event_counts = get_theta2_from_dl3(
                        dl3_path, 
                        target_coords=target_coords, 
                        theta2_axis=theta2_axis, 
                        theta_cut=theta_cut
                        )
                    if tt == 21:
                        ax13[0].set_title('tel1')
                        plot_theta2_dl3(
                            ax=ax13[0], 
                            theta2_axis=theta2_axis, 
                            counts_on=counts_on, 
                            counts_off=counts_off, 
                            alpha=alpha, 
                            theta_cut=theta_cut, 
                            event_counts=event_counts
                            )
                    elif tt == 22:
                        ax13[1].set_title('tel2')
                        plot_theta2_dl3(
                            ax=ax13[1], 
                            theta2_axis=theta2_axis, 
                            counts_on=counts_on, 
                            counts_off=counts_off, 
                            alpha=alpha, 
                            theta_cut=theta_cut, 
                            event_counts=event_counts
                            )                
                    else:
                        ax13[2].set_title('stereo')
                        plot_theta2_dl3(
                            ax=ax13[2], 
                            theta2_axis=theta2_axis, 
                            counts_on=counts_on, 
                            counts_off=counts_off, 
                            alpha=alpha, 
                            theta_cut=theta_cut, 
                            event_counts=event_counts
                            )
            if is_dist:
                histograms, histograms_diff, _, _, livetimes, survived_ped, bins = load_distributions_sst1m(dist_path=dist_path, dl3_path=dl3_path)

            # PLOTTING

            # Distributions
            ax4[0].set_title('Differential rates of DL1 events (per wobble)')
            ax4[0].set_ylabel('Differential rate')
            ax4[0].set_xlabel('Hillas Intensity [p.e.]')
            ax4[1].set_xlabel('livetime [s]')
            ax4[1].set_title('Livetime of each wobble')

            ax5[0].set_title('Differential rates of DL1 events (per wobble)')
            ax5[0].set_ylabel('Differential rate')
            ax5[0].set_xlabel('Hillas Intensity [p.e.]')
            ax5[1].set_xlabel('livetime [s]')
            ax5[1].set_title('Livetime of each wobble')

            ax6[0].set_title('Differential rates of DL1 events (per wobble)')
            ax6[0].set_ylabel('Differential rate')
            ax6[0].set_xlabel('Hillas Intensity [p.e.]')
            ax6[1].set_xlabel('livetime [s]')
            ax6[1].set_title('Livetime of each wobble')

            if is_dist and (tt==21):
                ax4[0].plot(bins[:-1], histograms_diff.T, alpha=0.5)
                ax4[0].set_xscale('log')
                ax4[0].set_yscale('log')
                ax4[0].axvline(50)
                ax4[0].set_xlim([10, 10**5])
                ax4[0].grid()

                ax4[1].hist(livetimes, bins=20)
                ax4[1].grid()
                fig4.suptitle(tel, fontsize=16)

            if is_dist and (tt==22):
                ax5[0].plot(bins[:-1], histograms_diff.T, alpha=0.5)
                ax5[0].set_xscale('log')
                ax5[0].set_yscale('log')
                ax5[0].axvline(50)
                ax5[0].set_xlim([10, 10**5])
                ax5[0].grid()

                ax5[1].hist(livetimes, bins=20)
                ax5[1].grid()
                fig5.suptitle(tel, fontsize=16)

            if is_dist and (tt==0):
                ax6[0].plot(bins[:-1], histograms_diff.T, alpha=0.5)
                ax6[0].set_xscale('log')
                ax6[0].set_yscale('log')
                ax6[0].axvline(50)
                ax6[0].set_xlim([10, 10**5])
                ax6[0].grid()

                ax6[1].hist(livetimes, bins=20)
                ax6[1].grid()
                fig6.suptitle(tel, fontsize=16)

            # Zenith angles
            fig7.suptitle("Distribution of zenith angles", fontsize=16)
            ax7[0].set_xlabel('zenith [deg]')
            ax7[0].set_ylabel('N reconstructed events')
            ax7[1].set_xlabel('local_time [s]')
            ax7[1].set_ylabel('zenith [deg]')

            if is_dl2 and (tel != "stereo"):
                ax7[0].grid()
                ax7[0].grid()
                ax7[0].legend()
                ax7[1].grid()
                ax7[1].invert_yaxis()
                ax7[1].legend()

            # CoGs
            # DL1 mono, survived cleaning
            fig8.suptitle("CoG of DL1 mono events", fontsize=16)
            if is_dl1 and (tt == 21):
                ax8[0].set_title('tel_021')
            if is_dl1 and (tt == 22):
                ax8[1].set_title('tel_022')
            
            # DL2 mono
            # Gammas
            fig9.suptitle("CoG of DL2 mono events, gammaness > 0.7", fontsize=16)
            if is_dl2 and (tt == 21):
                ax9[0].set_title(tel)
            if is_dl2 and (tt == 22):
                ax9[1].set_title(tel)

            # DL1 stereo, survived cleaning
            fig10.suptitle("CoG of DL1 stereo events", fontsize=16)
            fig12.suptitle("Time differences between stereo events", fontsize=16)
            ax10[0].set_title('tel1')
            ax10[1].set_title('tel2')
            ax12[0].set_xlabel('DT [ns]')
            ax12[1].set_xlabel('local time [s]')
            ax12[1].set_ylabel('DT [ns]')
            if is_dl1 and (tt == 0):
                # WR dt
                ax12[0].set_xlim([-1200, 1200])
                ax12[0].grid()
                ax12[0].set_yscale('log')
                ax12[1].set_ylim([-1200, 1200])

            # Moon
            ax11.set_title('Moon')  
            ax11.set_ylabel('deg')
            ax11.set_xlabel('Time [UTC]')
            if is_dl2:
                new_times = []
                for i in range(len(ax11.get_xticklabels())):
                    time_unix = ax11.get_xticklabels()[i].get_position()[0]
                    new_times.append(Time(time_unix, format='unix', scale='utc').isot.split('T')[1].split('.')[0])
                ax11.set_xticks(ax11.get_xticks())
                ax11.set_xticklabels(new_times) #, rotation='vertical')
                ax11.grid()
                ax11.legend()
                
            # changing tick labels to UTC time
            ax.set_title('DL1 events (survived cleaning)', fontsize=16)
            ax.set_xlabel('Time [UTC]')
            ax.set_ylabel('Trigger rate [Hz]')
            if is_dl1:
                new_times = []
                for i in range(len(ax.get_xticklabels())):
                    time_unix = ax.get_xticklabels()[i].get_position()[0]
                    new_times.append(Time(time_unix, format='unix', scale='utc').isot.split('T')[1].split('.')[0])
                ax.set_xticks(ax.get_xticks())
                ax.set_xticklabels(new_times) #, rotation='vertical')
                ax.grid()
                ax.legend()

            # changing tick labels to UTC time
            ax1.set_title('DL1 events (zoom)', fontsize=16)
            ax1.set_xlabel('Time [UTC]')
            ax1.set_ylabel('Trigger rate [Hz]')
            if is_dl1:
                new_times = []
                for i in range(len(ax1.get_xticklabels())):
                    time_unix = ax1.get_xticklabels()[i].get_position()[0]
                    new_times.append(Time(time_unix, format='unix', scale='utc').isot.split('T')[1].split('.')[0])
                ax1.set_xticks(ax1.get_xticks())
                ax1.set_xticklabels(new_times) #, rotation='vertical')
                ax1.grid()
                ax1.legend()
                ax1.set_ylim([0, 2*max(median1)])
                

            # changing tick labels to UTC time
            ax2.set_title('DL2 events', fontsize=16)
            ax2.set_xlabel('Time [UTC]')
            ax2.set_ylabel('Event rate [Hz]')
            if is_dl2:
                new_times = []
                for i in range(len(ax2.get_xticklabels())):
                    time_unix = ax2.get_xticklabels()[i].get_position()[0]
                    new_times.append(Time(time_unix, format='unix', scale='utc').isot.split('T')[1].split('.')[0])
                ax2.set_xticks(ax2.get_xticks())
                ax2.set_xticklabels(new_times) #, rotation='vertical')
                ax2.set_ylim([0, 2*max(median2)])
                ax2.grid()
                ax2.legend()

            ax3.set_xlabel('Time [UTC]')
            ax3.set_ylabel('NSB rate [MHz]')
            ax3.set_title("Night Sky Background", fontsize=16)
            if is_dl1:
                ax3.grid()
                ax3.legend()
    
        logging.info('All done, saving figures in temporary png files.')
        fig.savefig(outpath+'/trigger_rates.png', dpi=250)
        fig1.savefig(outpath+'/trigger_rates_zoom.png', dpi=250)
        fig2.savefig(outpath+'/reco_rates_zoom.png', dpi=250)
        fig3.savefig(outpath+'/nsb.png', dpi=250)
        fig4.savefig(outpath+'/diff_rates_1.png', dpi=250)
        fig5.savefig(outpath+'/diff_rates_2.png', dpi=250)
        fig6.savefig(outpath+'/diff_rates_stereo.png', dpi=250)
        fig7.savefig(outpath+'/zenith_angles.png', dpi=250)
        fig8.savefig(outpath+'/cog_dl1_mono.png', dpi=250)
        fig9.savefig(outpath+'/cog_dl2_mono_photons.png', dpi=250)
        fig10.savefig(outpath+'/cog_dl1_stereo.png', dpi=250)
        fig11.savefig(outpath+'/moon.png', dpi=250)
        fig12.savefig(outpath+'/stereo_dt.png', dpi=250)
        fig13.savefig(outpath+'/theta2.png', dpi=250)
        plt.close("all")

        # store pdf
        logging.info('Combining png images in the final pdf file.')
        img1 = prepare_to_pdf(outpath+'/trigger_rates.png')
        img2 = prepare_to_pdf(outpath+'/trigger_rates_zoom.png')
        img3 = prepare_to_pdf(outpath+'/reco_rates_zoom.png')
        img4 = prepare_to_pdf(outpath+'/nsb.png')
        img5 = prepare_to_pdf(outpath+'/diff_rates_1.png')
        img6 = prepare_to_pdf(outpath+'/diff_rates_2.png')
        img7 = prepare_to_pdf(outpath+'/diff_rates_stereo.png')
        img8 = prepare_to_pdf(outpath+'/zenith_angles.png')
        img9 = prepare_to_pdf(outpath+'/cog_dl1_mono.png')
        img10 = prepare_to_pdf(outpath+'/cog_dl2_mono_photons.png')
        img11 = prepare_to_pdf(outpath+'/cog_dl1_stereo.png')
        img12 = prepare_to_pdf(outpath+'/moon.png')
        img13 = prepare_to_pdf(outpath+'/stereo_dt.png')
        img14 = prepare_to_pdf(outpath+'/theta2.png')
        image_list = [img2, img3, img4, img5, img6, img7, img8, img9, img10, img11, img13, img12, img14]

        pdf_file = 'night_summary_' + str(date) + '_' +source+'.pdf'
        img1.save(outpath + '/' + pdf_file, "PDF", save_all=True, resolution=100.0, append_images=image_list)

        # remove all temporary png images
        logging.info('Removing all temporary png images.')
        for f in glob.glob(outpath+'/*.png'):
            os.remove(f)

if __name__ == "__main__":
    main()