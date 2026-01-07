#!/usr/bin/env python

"""
Script to calculate observability of any source from the catalog
for a single particular night. If one wants to quickly see the
observability of a single source, which does not need to be put
in the catalog (some transient for example), the coordinates may 
be specified independently. If no date is specified, the following
night is the default. Any source in the catalog can be removed from 
plotting/calculation if obsplan_show: false. The script plots source
visibility and moon distance and stores the figure in the output dir.
It also prints observability windows with moon distance and altitude
taken into account. All times are in local time (CEST/CET).

Usage:

$> python sst1mpipe_obsplan.py
--date 20240810
--source-catalog sst1m_source_catalog.json
--output-dir ./
--config sst1mpipe_config.json
--min-moon-dist 50.
--max-zenith 50.
--ra 350.
--dec 75.
"""

import argparse
import logging
import sys
from datetime import datetime

import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz, SkyCoord, get_moon, get_sun
from astropy.time import Time

import sst1mpipe
from sst1mpipe.io import check_outdir, load_config, load_source_catalog
from sst1mpipe.utils import get_location, get_moon_phase


def parse_args():
    parser = argparse.ArgumentParser(description="Plots observability of selected catalog sources for a single night of interest and calculates the observability windows.")

    parser.add_argument(
                        '--date',
                        dest='date',
                        default=None,
                        type=str,
                        help='Date of the observing night (YYYYMMDD). If not specified coming night is default (but be careful when running it after midnight!).',
                        )
    parser.add_argument(
                        '--out-dir', '-o',
                        dest='out_dir',
                        type=str,
                        help='Output directory to store figure with path of sources in the sky.',
                        default='./'
                        )
    parser.add_argument('--source-catalog', action='store', type=str,
                        dest='source_cat',
                        help='Path to a json file with catalog of all observed SST1M targets.',
                        )

    parser.add_argument('--config', '-c', action='store', type=str,
                        dest='config_file',
                        help='Path to a configuration file (to get the telescope location).',
                        )
    parser.add_argument(
                        '--min-moon-dist',
                        dest='min_moon',
                        default=50.,
                        type=float,
                        help='Minimum Moon distance (deg) for calculation of observability window.',
                        )
    parser.add_argument(
                        '--max-zenith',
                        dest='max_zenith',
                        default=50.,
                        type=float,
                        help='Maximum zenith angle (deg) for calculation of observability window.',
                        )
    parser.add_argument(
                        '--ra',
                        dest='ra',
                        type=float,
                        help='Right ascention (deg) for a quick check of a source of interest if one does not want to have it in catalog.',
                        )                    
    parser.add_argument(
                        '--dec',
                        dest='dec',
                        type=float,
                        help='Declination (deg) for a quick check of a source of interest if one does not want to have it in catalog.',
                        )   
    args = parser.parse_args()
    return args


def main():

    args = parse_args()
    date = args.date
    source_catalog_file = args.source_cat
    config_file = args.config_file
    min_moon_dist = args.min_moon * u.deg
    max_zenith = args.max_zenith * u.deg
    min_alt = 90.*u.deg - max_zenith
    out_dir = args.out_dir
    ra = args.ra
    dec = args.dec

    check_outdir(out_dir)

    logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers= [
                logging.FileHandler(out_dir+'/sst1mpipe_obsplan_'+str(sst1mpipe.__version__) + '.log', 'w+'),
                logging.StreamHandler(stream=sys.stdout)
                ]
        )
    logging.info('sst1mpipe version: %s', sst1mpipe.__version__)

    try:
        source_catalog = load_source_catalog(source_catalog_file)
    except Exception:
        logging.warning('Source catalog file not found!')
        source_catalog = {}
    if (ra is not None) and (dec is not None):
        if len(source_catalog) == 0:
            source_catalog = {'NEW': {'ra': 0}}
        source_catalog['NEW']['ra'] = ra
        source_catalog['NEW']['dec'] = dec
        source_catalog['NEW']['frame'] = 'icrs'
        source_catalog['NEW']['obsplan_show'] = True
    if len(source_catalog) == 0:
        logging.error('No source specified! Either provide a catalog file, or ra/dec of some source, or both.')

    if date is None:
        tonight_datetime = Time.now()
        tonight = tonight_datetime.datetime.strftime("%Y-%m-%d")
    else:
        tonight = str(date[:4]) + '-' + str(date[4:6]) + '-' + str(date[6:])
    print('\nCalculating source visibilities for night starting on', tonight, '\n')

    utc_shift = (datetime.now() - datetime.utcnow()).seconds / 3600 * u.hour
    if utc_shift > 1.5 * u.hour:
        scale = '[CEST]'
    else:
        scale = '[CET]'

    config = load_config(config_file)
    location = get_location(config=config, tel='tel_021')
    midnight = Time(tonight + ' 00:00:00') + 1 * u.day
    time_span = midnight + np.linspace(-8, 8, 100) * u.hour
    time_span_sparse = midnight + np.linspace(-8, 8, 10) * u.hour
    frame_night = AltAz(obstime=time_span, location=location)
    frame_night_sparse = AltAz(obstime=time_span_sparse, location=location)

    sun_altaz = get_sun(time_span).transform_to(frame_night)
    moon_altaz = get_moon(time_span).transform_to(frame_night)
    phase_angle_moon = get_moon_phase(times=midnight, loc=location)

    mask_sunset = sun_altaz.alt < -0*u.deg
    mask_nautical = sun_altaz.alt < -12*u.deg
    mask_astronomical = sun_altaz.alt < -18*u.deg

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    plt.plot((time_span + utc_shift).datetime, moon_altaz.alt, color='k', ls='--', label='Moon (phase = '+str(round(phase_angle_moon.to_value(u.deg), 1))+' deg)')
    plt.fill_between((time_span + utc_shift).datetime, 0, 90, mask_sunset, color='0.5', zorder=0, alpha=0.5) # civil twilight
    plt.fill_between((time_span + utc_shift).datetime, 0, 90, mask_nautical, color='0.2', zorder=0, alpha=0.5) # nautical twilight
    plt.fill_between((time_span + utc_shift).datetime, 0, 90, mask_astronomical, color='k', zorder=0, alpha=0.5) # astronomical twilight

    time_span.format = 'iso'
    sunset = min(time_span[mask_sunset] + utc_shift)
    print(('\nSunset'.ljust(23) + scale + ': ').ljust(31), get_time_formated(sunset))

    if sum(mask_nautical):
        nautical_twilight_evening = min(time_span[mask_nautical] + utc_shift)
        print(('Nautical evening'.ljust(22) + scale + ': ').ljust(30), get_time_formated(nautical_twilight_evening))
    else:
        print(('Nautical evening'.ljust(22) + scale + ': ').ljust(30), 'None')

    if sum(mask_astronomical):
        astronomical_twilight_evening = min(time_span[mask_astronomical] + utc_shift)
        print(('Astronomical evening'.ljust(22) + scale + ': ').ljust(30), get_time_formated(astronomical_twilight_evening))
    else:
        print(('Astronomical evening'.ljust(22) + scale + ': ').ljust(30), 'None')

    if sum(mask_astronomical):
        astronomical_twilight_morning = max(time_span[mask_astronomical] + utc_shift)
        print(('\nAstronomical morning'.ljust(23) + scale + ': ').ljust(31), get_time_formated(astronomical_twilight_morning))
    else:
        print(('\nAstronomical morning'.ljust(23) + scale + ': ').ljust(31), 'None')

    if sum(mask_nautical):
        nautical_twilight_morning = max(time_span[mask_nautical] + utc_shift)
        print(('Nautical morning'.ljust(22) + scale + ': ').ljust(30), get_time_formated(nautical_twilight_morning))
    else:
        print(('Nautical morning'.ljust(22) + scale + ': ').ljust(30), 'None')

    sunrise = max(time_span[mask_sunset] + utc_shift)
    print(('Sunrise'.ljust(22) + scale + ': ').ljust(30), get_time_formated(sunrise))

    print('\nObservability windows for different sources. Min altitude ('+str(min_alt.to_value(u.deg))+' deg), min Moon distance ('+str(min_moon_dist.to_value(u.deg))+' deg) taken into account.')
    print('ASTRO .. obs during astronomical night, NAUT .. obs during nautical night.')
    for source in source_catalog:
        if source_catalog[source]['obsplan_show']:
            if source_catalog[source]['frame'] == 'icrs':
                target_coords = SkyCoord(
                    ra=source_catalog[source]['ra']*u.deg, 
                    dec=source_catalog[source]['dec']*u.deg, 
                    frame='icrs'
                    )
            elif source_catalog[source]['frame'] == 'galactic':
                target_coords = SkyCoord(
                    l=source_catalog[source]['l']*u.deg, 
                    b=source_catalog[source]['b']*u.deg, 
                    frame='galactic'
                    )
            source_altaz = target_coords.transform_to(frame_night)
            p = plt.plot((time_span + utc_shift).datetime, source_altaz.alt, label=source)

            # Moon distance
            source_altaz_sparse = target_coords.transform_to(frame_night_sparse)
            moon_distance = source_altaz_sparse.separation(get_moon(time_span_sparse).transform_to(frame_night_sparse))
            color = p[0].get_color()

            for i in range(len(time_span_sparse.datetime)):
                if source_altaz_sparse.alt[i] > 0:

                    if ((time_span_sparse[i]+0.5*u.hour + utc_shift) < max(time_span_sparse)) & ((time_span_sparse[i]+0.5*u.hour + utc_shift) > min(time_span_sparse)):
                        plt.text((time_span_sparse+0.5*u.hour + utc_shift).datetime[i], 
                        source_altaz_sparse.alt.to_value(u.deg)[i], 
                        str(int(moon_distance.to_value(u.deg)[i])), 
                        color=color
                        )

            # printing observability windows
            print('\n## '+source)
            moon_distance_detailed = source_altaz.separation(get_moon(time_span).transform_to(frame_night))
            mask_moon = moon_distance_detailed > min_moon_dist
            mask_zenith = source_altaz.alt > min_alt
            mask_astro = mask_moon * mask_zenith * mask_astronomical
            mask_naut = mask_moon * mask_zenith * mask_nautical

            if sum(mask_astro):
                print(  'ASTRO'.ljust(5), 
                        scale.ljust(10), 
                        get_time_formated(min(time_span[mask_astro] + utc_shift)), 
                        ' - ',
                        get_time_formated(max(time_span[mask_astro] + utc_shift))
                        )
            else:
                print('ASTRO'.ljust(5), 'None')
            if sum(mask_naut):
                print(  'NAUT'.ljust(5), 
                        scale.ljust(10), 
                        get_time_formated(min(time_span[mask_naut] + utc_shift)), 
                        ' - ',
                        get_time_formated(max(time_span[mask_naut] + utc_shift))
                        )
            else:
                print('NAUT'.ljust(5), 'None')

    #plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=45)
    myFmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.title(tonight)
    plt.ylim([0, 90])
    plt.xlim([min((time_span + utc_shift).datetime), max((time_span + utc_shift).datetime)])
    plt.xlabel('Time ' + scale)
    plt.ylabel('Altitude [deg]')
    plt.tight_layout()
    plt.grid()
    plt.legend()
    fig.savefig(out_dir + '/' + 'sst1m_observability_'+ tonight +'.png', dpi=250)
    plt.show()

def get_time_formated(time, format="hhmm"):

    time_hhmmss = time.to_string().split(' ')[1]
    time_hhmm = time_hhmmss.split(":")[0] + ":" + time_hhmmss.split(":")[1]
    if format == "hhmm":
        return time_hhmm
    elif format == "hhmmss":
        return time_hhmmss

if __name__ == "__main__":
    main()