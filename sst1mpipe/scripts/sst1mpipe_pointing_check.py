#!/usr/bin/env python

"""
A script to manually check the telescope pointing.
- Input is a single raw .fits.fz data file (containing single telescope data).
- Output is a png image showing a mean waveform std in each pixel together with 
position of brigth stars in the FoV

Usage:

$> python sst1mpipe_pointing_check.py
--input-file SST1M1_20240304_0012.fits.fz
--output-dir ./
--config sst1mpipe_config.json
--max-magnitude 8
--pointing-ra 85.0
--pointing-dec 25.0
--force-pointing

"""

import sst1mpipe
from sst1mpipe.io import (
    load_config,
)

from ctapipe.instrument import SubarrayDescription

from sst1mpipe.io.sst1m_event_source import SST1MEventSource

from sst1mpipe.utils import get_target

from sst1mpipe.utils import (
    simbad_query,
    get_location,
    get_tel_string,
)

import os
import sys
import argparse
import numpy as np
import logging

from os import path

import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time

from ctapipe.coordinates import CameraFrame
from ctapipe.visualization import CameraDisplay

import matplotlib.pyplot as plt


def parse_args():

    parser = argparse.ArgumentParser(description="Pointing check")

    # Required arguments
    parser.add_argument(
                        '--input-file', '-f', type=str,
                        dest='input_file',
                        help='Path to the simtelarray or data file',
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
                        help='Path to store the output DL1 file',
                        default='./'
                        )

    parser.add_argument(
                        '--max-magnitude', '-m', type=str,
                        dest='max_mag',
                        help='Maximum B magnitude of a star in Simbad conesearch.',
                        default=8
                        )

    # Temporary
    parser.add_argument(
                        '--pointing-ra', '-r', type=float,
                        dest='ra',
                        help='Pointing RA (deg)',
                        default=None
                        )

    parser.add_argument(
                        '--pointing-dec', '-d', type=float,
                        dest='dec',
                        help='Pointing DEC (deg)',
                        default=None
                        )

    parser.add_argument(
                        '--force-pointing',
                        action='store_true',
                        help='Use pointing coordinates provided manualy by user even if there is a pointing info in the fits file.',
                        dest='force_pointing'
                        )

    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    input_file = args.input_file
    outdir = args.outdir
    ra = args.ra
    dec = args.dec
    max_mag = args.max_mag
    force_pointing = args.force_pointing

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers= [
            logging.FileHandler(outdir+'/sst1mpipe_pointing_check.log', 'w+'),
            logging.StreamHandler(stream=sys.stdout)
            ]
    )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    # Reading target name and assumed pointing ra,dec from the TARGET field
    # of the Events fits header
    target, ra_fits, dec_fits, wobble_fits = get_target(input_file)
    if (ra_fits is not None) & (dec_fits is not None):
        logging.info('Pointing info from the fits file: TARGET: ' + target + ', COORDS: ' + str(ra_fits) + ' ' + str(dec_fits) + ', WOBBLE: ' + wobble_fits)
        if force_pointing & (ra is not None) & (dec is not None):
            logging.warning('Using pointing info from manual input anyway (forced by user).')
            wobble = input_file.split('/')[-1].split('.')[0].split('_')[-1]
        else:
            ra = ra_fits
            dec = dec_fits
            wobble = wobble_fits
    else:
        logging.warning('Using pointing info from manual input.')
        wobble = input_file.split('/')[-1].split('.')[0].split('_')[-1]

    result_table = simbad_query(ra=ra, dec=dec, max_mag=max_mag)

    config = load_config(args.config_file)

    wobble_coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

    source = SST1MEventSource([input_file], max_events=None)
    
    ### This sould probably be moved to sst1m_event_source
    subarray_file = pkg_resources.resource_filename(
                        'sst1mpipe',
                        path.join(
                            'data',
                            'sst1m_array.h5'
                        )
                    )
    subarray = SubarrayDescription.from_hdf(subarray_file, focal_length_choice="EQUIVALENT")
    source._subarray = subarray

    std_samples = []
    time_all = []

    logging.info('Reading raw waveforms from the input file...')

    for i, event in enumerate(source):
        # NOTE: This needs to be changed in the future when event source hopefuly provides events with both telescope data
        if i == 0:
            tel = event.sst1m.r0.tels_with_data[0]

        r0data = event.sst1m.r0.tel[tel]
        event_type = event.sst1m.r0.tel[tel]._camera_event_type.value
        
        if event_type == 8:
            std_samples.append(np.std(r0data.adc_samples, axis=1))
            time = event.sst1m.r0.tel[tel].local_camera_clock/10**9
            time_all.append(time)

    std_samples = np.array(std_samples)
    time_all = np.array(time_all)

    tel_string = get_tel_string(tel, mc=False)
    location = get_location(config=config, tel=tel_string)

    index_mid = int(std_samples.shape[0]/2)
    obstime = Time(time_all, format='unix')
    horizon_frame = AltAz(location=location, obstime=obstime[index_mid])
    tel_pointing = wobble_coords.transform_to(horizon_frame)

    # Length of one file
    file_time = round((obstime[-1]-obstime[0]).to_value(u.s), 1)

    focal = float(config['telescope_equivalent_focal_length'][tel_string]) * u.m
    camera_frame = CameraFrame(
        focal_length=focal, 
        telescope_pointing=tel_pointing, 
        obstime=obstime[index_mid],
        location=location)

    std_samples_mean = np.mean(std_samples, axis=0)

    shift_x = 0
    shift_y = 0
    
    geometry = subarray.tel[tel].camera.geometry

    f,ax = plt.subplots(1, 1, figsize=(15,15))
    disp = CameraDisplay(geometry,ax=ax)
    disp.add_colorbar(ax=ax)
    disp.autoscale = False
    disp.norm.vmax = np.mean(std_samples_mean) + 1
    disp.norm.vmin = np.mean(std_samples_mean) - 2
    disp.image = std_samples_mean

    for star in result_table:
        star1 = SkyCoord(star['RA'] + ' ' + star['DEC'], frame='icrs', unit=(u.hourangle, u.deg))
        star1_altaz = SkyCoord(alt=star1.transform_to(horizon_frame).alt,
                            az=star1.transform_to(horizon_frame).az, 
                            frame=horizon_frame)
        star1_camera_frame = star1_altaz.transform_to(camera_frame)
        ax.plot(star1_camera_frame.x.to('m').value+shift_x, star1_camera_frame.y.to('m').value+shift_y,'*',color='black',markersize=9,label='star1')

    date = input_file.split('/')[-1].split('_')[1]
    obs_id = input_file.split('/')[-1].split('.')[0].split('_')[2]
    plt.title(str(date) + '_' + str(obs_id) + ' Time:' + str(file_time))

    outfile = 'pedestal_std_' + date + '_' + obs_id + '_' + wobble + '.png'
    plt.savefig(outdir + '/' + outfile)

    plt.close()

if __name__ == '__main__':
    main()
