"""
Copyright (c) 2024 SST-1M collaboration
Licensed under the 3-clause BSD style license.
"""

from ctapipe.io import read_table
from astropy.coordinates import AltAz, SkyCoord, EarthLocation, get_moon, get_sun, get_body
from ctapipe.coordinates import CameraFrame
from ctapipe.image import camera_to_shower_coordinates
from ctapipe.instrument import SubarrayDescription

import ctaplot
from astropy.time import Time
from astropy.table import Table
from astropy.io import fits
import astropy.constants as c
import numpy as np
import pandas as pd
import astropy.units as u
import h5py
import logging
import tables
import os
from astroquery.simbad import Simbad
import pkg_resources
from os import path


def get_target(file, force_pointing=False):
    """
    Extracts the target information from the string 
    stored in the TARGET field of the input file 
    Fits header. The string is expected to have the 
    following format \"Target,RA[in deg],DEC[in deg]\".
    Target may contain wobble information. Files stored 
    during the transition between wobbles must be flagged 
    with \'Transition\' string.

    Parameters
    ----------
    file: string
        Path to the raw fits data file
    force_poiting: bool
        If True, Transition flag is ignored and the
        file is processed anyway

    Returns
    -------
    target: string
    ra: float
        RA in degress
    dec: float
        DEC in degress
    wobble: string

    """

    with fits.open(file) as hdul:

        try:
            header = hdul["Events"].header
            # Delimiter should be hopefuly either ',' or '_'
            # The string in the TARGET field is expected in the form: target[]wobble[]ra[]dec, 
            # but targetwobble[]ra[]dec, and also target_wobble[]ra[]dec should work as well
            pointing_string = header['TARGET']
            logging.info('TARGET field: ' + pointing_string)
            if pointing_string == 'Transition' and not force_pointing:
                logging.info('Transition to the next wobble, not stable pointing direction, FILE SKIPPED.')
                hdul.close()
                exit()
            if pointing_string.count('_') > 1:
                delimiter = '_'
            elif pointing_string.count(',') > 1:
                delimiter = ','
            else:
                logging.warning('Wrong format of coordinates in the fits header, unknown delimiter')
                target, ra, dec, wobble = None, None, None, None
                return target, ra, dec, wobble

            target = pointing_string.split(delimiter)[0]
            try:
                if len(pointing_string.split(delimiter)) == 4:
                    ra = float(pointing_string.split(delimiter)[2])
                    dec = float(pointing_string.split(delimiter)[3])
                elif len(pointing_string.split(delimiter)) == 3:
                    ra = float(pointing_string.split(delimiter)[1])
                    dec = float(pointing_string.split(delimiter)[2])
                else: 
                    logging.warning('Wrong format of coordinates in the fits header. Field with either 3 or 4 entries is expected.')
                    ra, dec = None, None
            except ValueError:
                logging.warning('Wrong format of coordinates in the fits header, cannot convert to float!')
                ra, dec = None, None
            if 'W1' in pointing_string:
                wobble = 'W1'
            elif 'W2' in pointing_string:
                wobble = 'W2'
            elif 'W3' in pointing_string:
                wobble = 'W3'
            elif 'W4' in pointing_string:
                wobble = 'W4'
            else: wobble = 'UNDEF'
        except KeyError:
            logging.warning('TARGET field is not in the fits header! Cannot read pointing RA, DEC. Are you sure that this is a valid file with science data?')
            target, ra, dec, wobble = None, None, None, None
    return target, ra, dec, wobble



def get_nsb_levels_rates(config):

    nsb_levels={}
    nsb_rates={}
    telescopes=list(config['mean_charge_to_nsb_rate'].keys())
    for tel in telescopes:
        nsb_levels[tel] = sorted([x['mean_charge_bin_low'] for x in config['mean_charge_to_nsb_rate'][tel]])
        nsb_rates[tel] = {x['mean_charge_bin_low']: x['nsb_rate'] for x in config['mean_charge_to_nsb_rate'][tel]}
    return nsb_levels, nsb_rates


def get_stereo_method(config):
    """
    Reads method of coincidet event matching 
    from the config file.

    Parameters
    ----------
    config: dict

    Returns
    -------
    string
        Selected method of stereo event matching

    """
    stereo_method = config['stereo']['event_matching_method']
    if stereo_method not in ['SlidingWindow', 'WhiteRabbitClosest', 'SWATEventIDs']:
        logging.error("Unknown event mathing method set in cfg file! Available options: SlidingWindow, WhiteRabbitClosest, SWATEventIDs.")
        exit(-1)
    else:
        return stereo_method


def get_wr_timestamp(data):
    """
    Reads WR timestamps stored in two columns in the DL1
    table in seconds, and converts them in nanoseconds
    with high numerical precission.

    Parameters
    ----------
    data: pandas.DataFrame
        DL1 table

    Returns
    -------
    numpy.int64
        Precise timestamp in nanoseconds

    """

    S_TO_NS = np.int64(1e9)
    t1 = np.array(data['time_wr_full_seconds']).astype(np.int64) * S_TO_NS
    t2 = np.array(data['time_wr_frac_seconds']).astype(np.float64) * S_TO_NS
    t = t1 + t2.astype(np.int64)
    return t


def get_tel_string(tel, mc=True):
    """
    Makes string with telescope designation understandable
    for the rest of the analysis chain.

    Parameters
    ----------
    tel: int
        Telescope number as in
        event.sst1m.r0.tels_with_data
    mc: bool

    Returns
    -------
    string

    """

    if mc:
        if '1' in str(tel):
            telescope = 'tel_001'
        else:
            telescope = 'tel_002'
    else:
        if '1' in str(tel):
            telescope = 'tel_021'
        else:
            telescope = 'tel_022'
    return telescope


def get_location(config=None, tel=None):
    """
    Creates astropy object EarthLocation for
    given telescope from parameters stored in
    the config file.

    Parameters
    ----------
    config: dict
    tel: string

    Returns
    -------
    astropy.coordinates.EarthLocation

    """

    longitude = config["telescope_coords"][tel]["long_deg"] * u.deg
    latitude = config["telescope_coords"][tel]["lat_deg"] * u.deg
    height = config["telescope_coords"][tel]["height_m"] * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, height)

    return location


def simbad_query(ra=None, dec=None, max_mag=8):
    """
    Query simbad for the brightest stars in the 
    FoV.

    Parameters
    ----------
    ra: float
        Pointing right ascention, deg
    dec: float
        Pointing declination, deg
    max_mag: float
        Max stellar B magnitude in query region

    Returns
    -------
    astropy.table.Table
        Table of coordinates and magnitudes of
        the stars in the FoV. The table is clipped
        on the first 10 entries.
    """

    logging.info('Querying Simbad for list of bright sources around requested position...')
    Simbad.reset_votable_fields()
    Simbad.add_votable_fields('flux(B)')
    if dec > 0:
        simbad_line = 'region(circle, ICRS, ' + str(ra) + ' +' + str(dec) + ', 5d) & Bmag<'+str(max_mag)
    else:
        simbad_line = 'region(circle, ICRS, ' + str(ra) + ' ' + str(dec) + ', 5d) & Bmag<'+str(max_mag)
    logging.info('Simbad entry: ' + simbad_line)
    result_table = Simbad.query_criteria(simbad_line)
    result_table.sort('FLUX_B')
    if len(result_table) >= 10:
        result_table = result_table[:10]
    logging.info('First 10 lines in the results table (stars to be plotted):')
    logging.info(result_table)

    return result_table


def event_hillas_add_units(event):
    """
    This returns units to Hillas parameters to 
    event which resulted from 
    ctapipe.io.HDF5EventSource 

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    Returns
    -------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer
        Event with units added to Hillas 
        parameter
    """

    for tel in event.trigger.tels_with_trigger:
        event.dl1.tel[tel].parameters.hillas.length = event.dl1.tel[tel].parameters.hillas.length * u.m
        event.dl1.tel[tel].parameters.hillas.length_uncertainty = event.dl1.tel[tel].parameters.hillas.length_uncertainty * u.m
        event.dl1.tel[tel].parameters.hillas.phi = event.dl1.tel[tel].parameters.hillas.phi * u.deg
        event.dl1.tel[tel].parameters.hillas.psi = event.dl1.tel[tel].parameters.hillas.psi * u.deg
        event.dl1.tel[tel].parameters.hillas.r = event.dl1.tel[tel].parameters.hillas.r * u.m
        event.dl1.tel[tel].parameters.hillas.width = event.dl1.tel[tel].parameters.hillas.width * u.m
        event.dl1.tel[tel].parameters.hillas.width_uncertainty = event.dl1.tel[tel].parameters.hillas.width_uncertainty * u.m
        event.dl1.tel[tel].parameters.hillas.x = event.dl1.tel[tel].parameters.hillas.x * u.m
        event.dl1.tel[tel].parameters.hillas.y = event.dl1.tel[tel].parameters.hillas.y * u.m
        event.dl1.tel[tel].parameters.timing.slope = event.dl1.tel[tel].parameters.timing.slope * 1/u.m
    
    return event


def add_trigger_time(event, telescope=None):
    """
    Stores local_camera_clock provided with event source
    in event.trigger with nanosecond precision, see
    https://github.com/cta-observatory/ctapipe_io_nectarcam/issues/24
    Note: If we read the data using ctapipe.io.read_table, 
    the numerical precision is lost anyway. This is the reason
    storing WR timestamps in two columns directly in the DL1
    table.

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer
    telescope: int

    Returns
    -------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    """

    localtime = event.sst1m.r0.tel[telescope].local_camera_clock.astype(np.uint64)
    # assuming local_camera_clock in gps (gps = tai - 19s), gps scale does not exist in astropy
    # tai = utc + 37 s (this is not constant in time and depend on leap seconds)
    #event.trigger.time = Time(localtime * u.s + 19 * u.s, format='unix', scale='tai') - 37 * u.s

    # assuming local_camera_clock in tai and conversion to utc
    #event.trigger.time = Time(localtime * u.s, format='unix', scale='tai') - 37 * u.s

    # assuming local_camera_clock in utc
    #event.trigger.time = Time(localtime, format='unix', scale='utc')

    # Time in event.trigger.time is stored in seconds, but should have ns precision, see
    # https://github.com/cta-observatory/ctapipe_io_nectarcam/issues/24
    # But if we read the data, using ctapipe.io.read_table, the numerical precision is lost anyway
    # We assume tai scale
    S_TO_NS = np.uint64(1e9)
    full_seconds = localtime // S_TO_NS
    fractional_seconds = (localtime % S_TO_NS) / S_TO_NS
    event.trigger.time = Time(full_seconds, fractional_seconds, format='unix_tai')
    event.trigger.tel[telescope].time = event.trigger.time

    return event


def add_event_id(event, filename=None, event_number=0):
    """
    Fills event container with event_id if missing.
    This invented event_id is based on date, run id
    and event number within the file. First two of
    those are extracted from the filename.

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer
    filename: string
    event_number: int

    Returns
    -------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    """

    date = filename.split('/')[-1].split('_')[1]
    obs_id = date + filename.split('/')[-1].split('_')[2]

    if event.sst1m.r0.event_id > 0:
        event_id = event.sst1m.r0.event_id
    else:
        event_id = str(int(filename.split('/')[-1].split('_')[2])) + str(event_number).zfill(6)
        logging.warning('Event IDs are not stored in raw data. Replacing with event ID based on date and event count.')

    event.index.event_id = int(event_id)
    event.index.obs_id = int(obs_id)

    return event


def add_pointing_to_events(
    event, ra=None, dec=None, 
    telescope=None, location=None):
    """
    Fills array pointing and telescope pointing 
    fields of the event container. Pointing 
    altitude ant azimuth are also calculated 
    from poiting ra,dec and trigger time.

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer
    ra: float
    dec: float
    telescope: int
    location: astropy.coordinates.EarthLocation

    Returns
    -------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    """

    wobble_coords = SkyCoord(ra=float(ra) * u.deg, dec=float(dec) * u.deg, frame='icrs')
    horizon_frame = AltAz(obstime=event.trigger.time, location=location)
    tel_pointing = wobble_coords.transform_to(horizon_frame)
    event.pointing.tel[telescope].azimuth  = tel_pointing.az.to('rad')
    event.pointing.tel[telescope].altitude = tel_pointing.alt.to('rad')
    event.pointing.array_azimuth  = tel_pointing.az.to('rad')
    event.pointing.array_altitude = tel_pointing.alt.to('rad')
    event.pointing.array_ra  = wobble_coords.ra.to('rad')
    event.pointing.array_dec = wobble_coords.dec.to('rad')
    return event


def check_mc(file):
    """
    Checks the input DL1 HDF file and
    returns True of False if it's 
    MC or data file, respectively.

    Parameters
    ----------
    file: string
        Input DL1 file.

    Returns
    -------
    bool

    """

    try:
        mc = read_table(file, "/configuration/simulation/run")
        logging.info('This is MC file.')
        return True
    except tables.exceptions.NoSuchNodeError:
        logging.info('This is REAL DATA file.')
        return False


def get_avg_pointing(dl1):

    tel_ze = 90.-dl1['true_alt_tel'].mean()
    tel_az = dl1['true_az_tel'].mean()

    alt_range = dl1['true_alt_tel'].max() - dl1['true_alt_tel'].min()
    az_range = dl1['true_az_tel'].max() - dl1['true_az_tel'].min()
    logging.info('Average pointing: zenith: %f, azimuth: %f', tel_ze, tel_az)
    logging.info('Range of zeniths: %f', alt_range)
    logging.info('Range of azimuths: %f', az_range)
    return tel_ze, tel_az


def get_closest_rf_model(
        dl1, models_dir=None, nsb_level=None, 
        tel=None, config=None):
    """
    RF models are trained per bin in zenith angle,
    azimuth and NSB level. They are expected to be
    stored in subdirectories with the following 
    naming format: \'zeXX_azXX_nsbXX\'. This provides
    the path to the closest RF models (in [NSB,ze,az]
    space).

    Parameters
    ----------
    dl1: pandas.DataFrame or astropy.table.Table
    models_dir: str
        Path to stored trained RFs, or path to general 
        production directory where subdirectories with 
        models in the following naming format are 
        expected: \'zeXX_azXX_nsbXX\''
    nsb_level: float
        Average NSB level, usualy stored in meanQ
        column of the DL1 table
    tel: string
    config: dict

    Returns
    -------
    string

    """

    dir_content = os.listdir(models_dir)
    if len([i for i in dir_content if '.sav' in i]) == 0:
        tel_ze, tel_az = get_avg_pointing(dl1)
        models_tab = []
        for dir_ze_az in dir_content:
            try:
                model_coords = dir_ze_az.split('_')
                models_tab.append(
                    [model_coords[0].split('ze')[1], 
                    model_coords[1].split('az')[1],
                    model_coords[2].split('nsb')[1],
                    dir_ze_az
                    ])
            except:
                logging.warning('%s does not follow the naming convention for RF model subdirectories, skipping..', os.path.join(models_dir, dir_ze_az))
        models_tab = np.array(models_tab)
        logging.info('%d RF nodes found in %s.',len(models_tab), models_dir)
        if len(models_tab) == 0:
            logging.error('Finding of the closest RFs failed. Subdirectories with RFs probably does not follow required naming format: \'zeXX_azXX_nsbXX\'')

        # Find closest zenith and NSB
        # TODO: closest node in general using azimuth
        # Factor of 666 below is an arbitrary "high" number to force NSB level - we reconstruct for known NSB bins so exact match is expected
        if (config is not None) and (nsb_level is not None) and (tel is not None):
            nsb_levels, nsb_rates = get_nsb_levels_rates(config)
            for meanQ_low in reversed(nsb_levels[tel]):
                if nsb_level >= meanQ_low:
                    break
            logging.info('Average charge from pedestal events in the file is %f.', nsb_level)
            idx_closest = np.argmin(abs(models_tab[:, 0].astype(np.float64) - tel_ze) + 666*abs(models_tab[:, 2].astype(np.float64) - nsb_rates[tel][meanQ_low]))
        else:
            nsb_rate_default=100
            idx_closest = np.argmin(abs(models_tab[:, 0].astype(np.float64) - tel_ze) + 666*abs(models_tab[:, 2].astype(np.float64) - nsb_rate_default))

        closest_model_dir = os.path.join(models_dir, models_tab[idx_closest, 3])
        logging.info('Closest RF node found: %s', closest_model_dir)
        return closest_model_dir
    else:
        logging.info('There are RFs found directly in the directory provided: %s These will be used for the reconstruction without any further selection.', models_dir)
        return models_dir


def check_same_shower_fraction(dl2, energy_bins):

    logging.info('SANITY CHECK: In each bin in reconstructed energy a fraction of reused triggered events is shown.')
    logging.info('E.g.: [0.5 0.25 0.25 0 0 0 0 0 0 ...] means that half of all showers in given energy bin was used once, 25%% of showers was used 2x, 25%% 3x, and that\'s it.')
    n_reuse_bins = np.arange(1, 31, 1)
    logging.info('N reuse bins: %s', list(n_reuse_bins))
    energy = dl2['reco_energy']

    # Same simulated showers can be identified by true energy, which might be a bit dangerous, if the numeric precision is not high enough
    # to catch all showers which differ only by a tiny fraction of E_True. An alternative is to use obs_id together with event_id, because
    # the events generated by the same shower have the same event_id BUT THE LAST TWO DIGITS. The numbering starts from 100 (first shower, 
    # first reuse), i.e. it works up to number_reuse = 100. What happens above that is a big mystery.
    event_id = dl2['event_id']
    event_id -= event_id % +100
    dl2['shower_id'] = event_id

    # Result returned is fraction of events in each energy bin, which are used more than once
    result = []
    logging.info('Same shower fraction based on obs_id, shower_id (WARNING: It works only if shower_reuse <= 100!).')
    for i in range(len(energy_bins)-1):
        mask = (energy > energy_bins[i]) & (energy < energy_bins[i + 1])
        n_same_events = dl2[mask].groupby(['obs_id', 'shower_id']).size()
        #n_same_events = dl2[mask].groupby('true_energy').size()
        if len(dl2[mask]) > 0:
            total_events = len(dl2[mask])
            fractions = ["%.2E" % elem for elem in list(np.histogram(n_same_events, bins=n_reuse_bins)[0] / total_events)]
            fraction_used_more = sum(np.histogram(n_same_events, bins=n_reuse_bins)[0][1:] / total_events)
            logging.info('E_R [%.2f, %.2f] TeV, total events: %d, frac used N>1: %.2E, fracs: %s', energy_bins[i].value, energy_bins[i + 1].value, total_events, fraction_used_more, fractions)
            result.append([energy_bins[i].value, energy_bins[i + 1].value, fraction_used_more])
        else:
            result.append([energy_bins[i].value, energy_bins[i + 1].value, np.nan])
    result = np.array(result)
    return result


def mc_correct_shower_reuse(mc_table, histograms):
    """
    If we use CSCAT>20 in CORSIKA the number stored in mc['shower_reuse']
    is still 20, even though the total number of shower stored in the simtel 
    file is correct (Vladimir checked), and also the total number of thrown 
    events stored in histograms seems to be fine (Jakub checked). Here we 
    use the total number of events stored in histograms to correct for 
    shower_reuse.

    Parameters
    ----------
    mc_table: astropy.table.Table
        Table extracted from DL1/DL2 MC
        file: /configuration/simulation/run
    histograms: astropy.table.Table
        Table extracted from DL1/DL2 MC
        file: /simulation/service/shower_distribution

    Returns
    ------- 
    astropy.table.Table

    """

    for simulated_e_min in mc_table.group_by('energy_range_min').groups.aggregate(sum)['energy_range_min']:

        mask = mc_table['energy_range_min'] == simulated_e_min

        # NOTE: just be careful if you decide to cut in histograms in the future
        n_true_simulated = histograms['n_entries'][mask].sum()
        n_simulated = (mc_table['n_showers'][mask] * mc_table['shower_reuse'][mask]).sum()
    
        if n_true_simulated != n_simulated:
            logging.warning('Sum of histograms for production with e_min: %f doesn\'t match with n_showers*shower_reuse stored in the mc tab!!', simulated_e_min)
            n_reuse_new = histograms['n_entries'][mask][0]  / mc_table['n_showers'][mask]
            logging.warning('CORRECTION OF SHOWER_REUSE APPLIED. New value: %f', n_reuse_new[0])
            mc_table['shower_reuse'][mask] = n_reuse_new

    return mc_table


def correct_number_simulated_showers(
        simulated_event_info, mc_table=None, 
        histogram_table=None, e_min_cut=None):
    """
    Corrects the number of totaly simulated showers 
    after the cuts on true energy in DL1 and DL2.

    One may combine different MC production with 
    different E mins and cut on the minimum energy
    when producing DL1 or DL2 files. But is such case
    the number of simulated events and simulated event
    histograms, wchi are stored in DL1/DL2 MC are not 
    correct resulting in wrong spectral reweighting.

    This fucntion corrects the number of simulated events
    and performs a simple sanity check.

    Parameters
    ----------
    simulated_event_info: astropy.table.QTable
    mc_table: astropy.table.Table
    histogram_table: astropy.table.Table
    e_min_cut: float
        Cut on E nergy applied in DL1/DL2 files. 
        Expected to be in TeV.

    Returns
    ------- 
    astropy.table.QTable

    """

    logging.info('Number of different MC productions combined: %d', len(mc_table.group_by('energy_range_min').groups.aggregate(sum)))
    N_showers_new = []

    for simulated_e_min in mc_table.group_by('energy_range_min').groups.aggregate(sum)['energy_range_min']:

        mask = mc_table['energy_range_min'] == simulated_e_min

        if simulated_e_min < e_min_cut:

            logging.info('Correcting simulated e min in %d files out of %d files.', sum(mask), len(mask))
            N_showers_old = (mc_table['n_showers'][mask] * mc_table['shower_reuse'][mask]).sum()
            exponent = mc_table['spectral_index'][0] + 1
            int_old = 1/exponent*mc_table['energy_range_max'][mask][0]**exponent - 1/exponent*mc_table['energy_range_min'][mask][0]**exponent
            int_new = 1/exponent*mc_table['energy_range_max'][mask][0]**exponent - 1/exponent*e_min_cut**exponent
            N_showers_new.append(int(N_showers_old * int_new / int_old))
            logging.info("N of simulated showers in original interval (%f, %f): %d", mc_table['energy_range_min'][mask][0], mc_table['energy_range_max'][mask][0], N_showers_old)
            logging.info("N of simulated showers in new interval (%f, %f): %d", e_min_cut, mc_table['energy_range_max'][mask][0], int(N_showers_old * int_new / int_old))

            # rough sanity check
            closest_bin = np.argmin(abs(histogram_table['bins_energy'][0][:-1] - e_min_cut))
            mask_bins = histogram_table['bins_energy'][0][:-1] >= histogram_table['bins_energy'][0][closest_bin]
            histogram_sum = sum(sum(histogram_table['histogram'][mask].sum(axis=2)[:, mask_bins]))
            logging.info('SANITY CHECK: Exact number of simulated events above %.5f TeV: %d', round(histogram_table['bins_energy'][0][closest_bin], 5), histogram_sum)
            difference = (histogram_sum - int(N_showers_old * int_new / int_old)) / histogram_sum
            logging.info('Difference [%%/100]: %f', difference)

        # We need to add also the rest of the files with given simulated e min, which have correct mc table
        else:
            logging.info('For %d files out of %d the energy cut is lower or equal than simulated, and the number of total simulated events dont have to be corrected.', sum(mask), len(mask))
            N_showers_new.append((mc_table['n_showers'][mask] * mc_table['shower_reuse'][mask]).sum())

    simulated_event_info.n_showers = sum(N_showers_new)
    logging.info('Total corrected number of simulated events: %d', simulated_event_info.n_showers)

    return simulated_event_info


def energy_min_cut(file, config=None):

    mc = read_table(file, "/configuration/simulation/run")
    particle_type = get_primary_type(file)

    if particle_type == 0:
        min_energy = config["analysis"]["gamma_min_simulated_energy_tev"] 
    elif particle_type == 101:
        min_energy = config["analysis"]["proton_min_simulated_energy_tev"]

    telescopes = get_telescopes(file)

    for tel in telescopes:

        params = read_table(file, "/dl1/event/telescope/parameters/" + tel)
        try:
            images_tab = True
            images = read_table(file, "/dl1/event/telescope/images/" + tel)
        except tables.exceptions.NoSuchNodeError:
            images_tab = False

        if (min_energy > mc['energy_range_min']).any():
            logging.info('Cutting out all events with true energies below %f TeV, which is above the minimum simulated energy %f TeV.', min_energy, mc['energy_range_min'][0])
            logging.warning("Note that beyond this point the mc table with number of all thrown events and the histograms are no longer correct. These will be corrected for automaticaly in performance evaluation.")

            mask = params['true_energy'] >= min_energy
            params['min_true_energy_cut'] = min_energy
            masked_params = params[mask]
            if images_tab:
                masked_images = images[mask]

            logging.info("Telescope %s: %d events removed.", tel, len(params)-len(masked_params))
        else:
            masked_params = params.copy()
            masked_params['min_true_energy_cut'] = mc['energy_range_min'][0]
            if images_tab:
                masked_images = images.copy()

        # Only combination of both overwrite and append works like expected, i.e. overwrite only telescope parameters and the rest 
        # of the DL1 file remains the same
        # NOTE: Unfortunately, we cannot store units with serialize_meta, because these are stored somehow weirdly as a new table and ctapipe merge tool
        # then cannot merge the files and raise error...
        masked_params.write(file, path='/dl1/event/telescope/parameters/'+tel, overwrite=True, append=True) #, serialize_meta=True)
        if images_tab:
            masked_images.write(file, path='/dl1/event/telescope/images/'+tel, overwrite=True, append=True) #, serialize_meta=True)


def remove_stereo(features):
    if "HillasReconstructor_tel_impact_distance" in features:
        features.remove("HillasReconstructor_tel_impact_distance")
    if "HillasReconstructor_h_max" in features:
        features.remove("HillasReconstructor_h_max")
    #if "swgo_reco_tilted_impact" in features:
    #    features.remove("swgo_reco_tilted_impact")
    return features


def get_event_sample(params, max_events=None):
    """
    Randomly draw subsamble of events.

    Parameters
    ----------
    params: pandas.DataFrame
        Input DL1/DL2 table
    max_events: int
        Sample size

    Returns
    -------
    pandas.DataFrame

    """

    if max_events is not None:
        max_events = int(max_events)
        if (len(params) < max_events):
            logging.info("Diffuse gamma training sample: \
                Requested number of training events cannot \
                be fulfilled, too few events in training DL1 files.")
        else:
            params = params.copy().sample(n=max_events)
        logging.info("Diffuse gamma training sample: N training gammas = %d", len(params))
    return params


def mix_gamma_proton(params_gammas, params_protons, max_events=None, gp_training_ratio=1.0):

    if max_events is not None:
        max_events = int(max_events)
        if (len(params_gammas) < gp_training_ratio/(gp_training_ratio + 1) * max_events):
            logging.info("GH classifier training: Requested number of training gammas cannot be fulfilled, too few events in training DL1 files.")
            N_gammas = len(params_gammas)
            if len(params_protons) < N_gammas / gp_training_ratio:
                logging.info("GH classifier training: Requested number of training protons cannot be fulfilled, too few events in training DL1 files.") 
                N_protons = len(params_protons)
            else:
                N_protons = int(N_gammas / gp_training_ratio)

        elif len(params_protons) < 1 / (1 + gp_training_ratio) * max_events:
            logging.info("GH classifier training: Requested number of training protons cannot be fulfilled, too few events in training DL1 files.")
            N_protons = len(params_protons)
            if len(params_gammas) < N_protons * gp_training_ratio:
                logging.info("GH classifier training: Requested number of training gammas cannot be fulfilled, too few events in training DL1 files.") 
                N_gammas = len(params_gammas)
            else:
                N_gammas = int(N_protons * gp_training_ratio)
        else:
            N_gammas = int(gp_training_ratio/(gp_training_ratio + 1) * max_events)
            N_protons = int(1 / (1 + gp_training_ratio) * max_events)

        params_gammas = params_gammas.copy().sample(n=N_gammas)
        params_protons = params_protons.copy().sample(n=N_protons)
        logging.info("GH classifier training: N training protons = %d, N training gammas = %d", N_protons, N_gammas)

    mixed = pd.concat([params_gammas, params_protons], ignore_index=True)
    mixed = mixed.sample(frac=1).reset_index(drop=True)
    return mixed


def correct_true_image(event):
    """
    For an unknown reason, event.simulation.tel[tel].true_image
    is sometimes None, which kills the rest of the script
    and simulation histogram is not saved. Here we repace 
    it with an array of zeros.

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    Returns
    -------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    """

    for tel in event.trigger.tels_with_trigger:
        if event.simulation.tel[tel].true_image is None:
            event.simulation.tel[tel].true_image = np.zeros(1296, dtype='int32')

    return event


def swap_modules_59_88(event, tel=None):
    """
    Swaps pixel R1 waveforms in two wrongly
    connected modules of tel2 after camera refurbishment
    in 2023. Pixel numbering is based on 
    test/resources/camera_config.cfg

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer
    tel: int

    Returns
    -------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    """

    # module 59
    mask59 = np.zeros(1296, dtype=bool)
    mask59[1029] = True
    mask59[1098:1102+1] = True
    mask59[1133:1134+1] = True
    mask59[1064:1067+1] = True
    waveform_59 = event.r1.tel[tel].waveform[mask59, :]
    
    # module 88
    mask88 = np.zeros(1296, dtype=bool)
    mask88[1103] = True
    mask88[1165:1169+1] = True
    mask88[1194:1195+1] = True
    mask88[1135:1138+1] = True
    waveform_88 = event.r1.tel[tel].waveform[mask88, :]
    
    event.r1.tel[tel].waveform[mask59] = waveform_88
    event.r1.tel[tel].waveform[mask88] = waveform_59   

    return event


def remove_bad_pixels(event, config=None):
    """
    Fills bad pixel waveforms with zeros and 
    flags them in proper containers. Charges in 
    these pixels are then interpolated using method 
    set in cfg: invalid_pixel_handler_type
    Default is NeighborAverage, but can be turned 
    off with 'null'

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer
    config: dict

    Returns
    -------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    """

    if "bad_pixels" in config["analysis"]:
        for tel in event.trigger.tels_with_trigger:
            tel_name = "tel_"+str(tel).zfill(3)
            if tel_name in config["analysis"]["bad_pixels"]:
                if len(config["analysis"]["bad_pixels"][tel_name]):

                    mask_bad = np.zeros(1296)
                    mask_bad[config["analysis"]["bad_pixels"][tel_name]] = 1
                    mask_bad = mask_bad.astype(bool)

                    N_samples = event.r0.tel[tel].waveform[0].shape[1]
                    event.r0.tel[tel].waveform[0][mask_bad] = np.zeros(N_samples)
                    event.r1.tel[tel].waveform[mask_bad] = np.zeros(N_samples)
                    event.simulation.tel[tel].true_image[mask_bad] = 0
                    event.mon.tel[tel].pixel_status['hardware_failing_pixels'] = np.array([mask_bad])
                    event.mon.tel[tel].pixel_status['flatfield_failing_pixels'] = np.array([mask_bad])
                    event.mon.tel[tel].pixel_status['pedestal_failing_pixels'] = np.array([mask_bad])

    return event


def check_output_dl1(file):
    """
    Checks if simulated shower distributions
    are properly stored in the DL1 MC file,
    which usualy means that the file is OK.

    Parameters
    ----------
    file: string

    Returns
    -------

    """

    try:
        hist = read_table(file, "/simulation/service/shower_distribution")
    except tables.exceptions.NoSuchNodeError:
        logging.error("EOFError reading simtel file. The file might be truncated and resulting DL1 file cannot be produced.")
        os.remove(file)


def add_features(data):
    """
    Add some extra parameters in the DL1 table. 

    Parameters
    ----------
    data: astropy.table.Table

    Returns
    -------
    astropy.table.Table 

    """

    data['log_camera_frame_hillas_intensity'] = np.log10(data['camera_frame_hillas_intensity'])
    data['camera_frame_hillas_wl'] = data['camera_frame_hillas_width']/data['camera_frame_hillas_length']

    return data


def add_log_true_energy(data):
    """
    Add log of the true energy in 
    the DL1 MC table. 

    Parameters
    ----------
    data: astropy.table.Table

    Returns
    -------
    astropy.table.Table 

    """

    data['log_true_energy'] = np.log10(data['true_energy'])
    return data


def add_pointing_mc(
        data, input_file=None, tel=None):
    """
    Add pointing altitude and azimuth in 
    DL1 table

    Parameters
    ----------
    data: astropy.table.Table
    input_file: string
    tel: string

    Returns
    -------
    astropy.table.Table 

    """

    string = "/dl1/monitoring/telescope/pointing/" + tel
    pointing = read_table(input_file, string)

    merged_table = data.copy()
    
    merged_table['true_az_tel'] = np.ones(len(merged_table)) * (pointing['azimuth'][0] * pointing['azimuth'].unit).to(u.deg)
    merged_table['true_alt_tel'] = np.ones(len(merged_table)) * (pointing['altitude'][0] * pointing['altitude'].unit).to(u.deg)
    
    return merged_table


def add_disp(paramtable, config=None, telescope=None):
    """
    Add DISP parameters in the DL1 MC table. 

    Parameters
    ----------
    paramtable: astropy.table.Table
    config: dict
    telescope: string

    Returns
    -------
    astropy.table.Table 

    """

    merged_table = paramtable
    
    horizon_frame = get_horizon_frame(config=config, telescope=telescope)

    true_coordinates = get_event_pos_in_camera(merged_table, horizon_frame, true_pos=True)
    merged_table['true_camera_x'] = true_coordinates[0]
    merged_table['true_camera_y'] = true_coordinates[1]
    
    # add disp parameters
    disp_parameters = add_disp_to_parameters(merged_table)
    
    merged_table['disp_dx'] = disp_parameters[0]
    merged_table['disp_dy'] = disp_parameters[1]
    merged_table['disp_norm'] = disp_parameters[2]
    merged_table['disp_angle'] = disp_parameters[3].value * u.rad
    merged_table['disp_sign'] = disp_parameters[4]
    
    return merged_table


def add_miss(paramtable):
    """
    Add miss parameter in the DL1 MC table. 

    Parameters
    ----------
    paramtable: astropy.table.Table

    Returns
    -------
    astropy.table.Table 

    """

    # NOTE: Be careful here, camera_frame_hillas_psi is in units of degrees
    disp, miss = camera_to_shower_coordinates(
        paramtable['true_camera_x'], 
        paramtable['true_camera_y'], 
        paramtable['camera_frame_hillas_x'], 
        paramtable['camera_frame_hillas_y'], 
        paramtable['camera_frame_hillas_psi'].value * u.deg
        )
    paramtable['miss'] = miss
    return paramtable


def add_true_impact(
        params, input_file=None, tel=None):
    """
    Add true impact parameter to DL1 MC table

    Parameters
    ----------
    params: astropy.table.Table
    input_file: string
    tel: string

    Returns
    -------
    astropy.table.Table 

    """

    array_info = read_table(input_file, "/configuration/instrument/subarray/layout")

    if tel == 'tel_001': index = 0
    else:  index = 1
    params['true_tel_impact_distance'] = np.sqrt((array_info['pos_x'][index] - params['true_core_x'])**2 + (array_info['pos_y'][index] - params['true_core_y'])**2) * u.m
    return params


def event_selection(data, config=None):
    """
    Performs event selection based on the cuts
    in the config file. It works on DL1/DL2,
    mono/stereo and MC/data

    Parameters
    ----------
    data: astropy.table.Table
    config: dict

    Returns
    -------
    astropy.table.Table 
        Table after quality cuts

    """

    events_selection = config["event_selection"]

    mask = np.ones(len(data), dtype=bool)

    logging.info('Application of selection cuts:')
    for key, (lower_limit, upper_limit) in events_selection.items():
        # This is because we use this very same function on event selection for DL2, and in stereo,
        # we have in DL2 these features with names changes so that we know which telescope they 
        # belong to.
        if key + "_tel1" in data.keys():
            for tel in ["tel1", "tel2"]:
                logging.info(f'{key + "_" + tel}: [{lower_limit}, {upper_limit}]')
                mask &= (data[key + "_" + tel] >= lower_limit) & (data[key + "_" + tel] <= upper_limit)
        elif key + "_tel21" in data.keys():
            for tel in ["tel21", "tel22"]:
                logging.info(f'{key + "_" + tel}: [{lower_limit}, {upper_limit}]')
                mask &= (data[key + "_" + tel] >= lower_limit) & (data[key + "_" + tel] <= upper_limit)
        else:
            logging.info(f'{key}: [{lower_limit}, {upper_limit}]')
            mask &= (data[key] >= lower_limit) & (data[key] <= upper_limit)

    return data[mask]


def get_finite(data, config=None, stereo=False):

    mask = np.ones(len(data), dtype=bool)

    N_selected = len(data[mask])
    features = set(config["energy_regression_features"] + 
            config["particle_classification_features"] + 
            config["disp_regression_features"] + 
            config["disp_classification_features"]
            )

    if not stereo:
        logging.info('Checking if feature values are finite excluding stereo params')
        features = remove_stereo(features)
    else:
        logging.info('Checking if feature values are finite INCLUDING STEREO PARAMS')

    for key in features:
        try:
            mask &= np.isfinite(data[key])
        except:
            logging.warning('{} column not in data.'.format(key))

    N_finite = len(data[mask])

    if N_selected > 0:
        if (N_finite/N_selected < 0.9) & (~stereo):
            logging.warning("There is more than 10% infinite values in features!")
    else:
        logging.warning("There is zero finite values in the file!")

    return data[mask]


def get_telescopes(input_file, data_level="dl1"):
    """
    Returns list of telescopes which have parameters
    table stored in the input DL1 or DL2 HDF file. 

    Parameters
    ----------
    input_file: string
        Input DL1 or DL2 file.
    data_level: string

    Returns
    -------
    list of strings
        List of telescopes in the file.

    """

    path=data_level+'/event/telescope/parameters'

    try:
        with h5py.File(input_file,'r') as f:
            telescopes = list(f[path].keys())
        return telescopes
    except KeyError:
        logging.error('No ' + path + ' table in the file: ' + input_file)
        return []


def get_primary_type(file):

    mc = read_table(file, '/simulation/event/subarray/shower')
    particle_type = mc['true_shower_primary_id'][0]
    return particle_type


def get_event_pos_in_camera(data, horizon_frame, true_pos=True):

    array_pointing = SkyCoord(
        alt=clip_alt(data['true_alt_tel']),
        az=data['true_az_tel'],
        frame=horizon_frame,
    )

    if true_pos:
        event_direction = SkyCoord(
            alt=clip_alt(data['true_alt']), az=data['true_az'], frame=horizon_frame
        )
    else:
        event_direction = SkyCoord(
            alt=clip_alt(data['reco_alt']), az=data['reco_az'], frame=horizon_frame
        )

    # NOTE: This should be checked. There is also effective_focal_length. Which one is better?
    focal = data['equivalent_focal_length']

    camera_frame = CameraFrame(focal_length=focal, telescope_pointing=array_pointing)

    camera_pos = event_direction.transform_to(camera_frame)
    return camera_pos.x, camera_pos.y


def camera_to_altaz(pos_x, pos_y, focal, pointing_alt, pointing_az, config=None, telescope=None, times=None):

    horizon_frame = get_horizon_frame(config=config, telescope=telescope, times=times)
    pointing_direction = SkyCoord(
        alt=clip_alt(pointing_alt), az=pointing_az, frame=horizon_frame
    )
    camera_frame = CameraFrame(
        focal_length=focal, telescope_pointing=pointing_direction
    )
    camera_coord = SkyCoord(pos_x, pos_y, frame=camera_frame)
    horizon = camera_coord.transform_to(horizon_frame)

    return horizon


def get_horizon_frame(config=None, telescope=None, times=None):

    location = get_location(config=config, tel=telescope)

    # This seems to have no effect
    if times is not None:
        obstime = times
    else:
        obstime = Time("2022-11-01T02:00")

    horizon_frame = AltAz(location=location, obstime=obstime)

    return horizon_frame


def get_survived_ped_fraction(dl1_files, logs=None, tel=None):

    zenith_angles = []
    times = []
    ped_fractions = []
    wobble = []
    filenames = []
    runs = []

    for file in dl1_files:
        try:
            zenith_angle = 90. - read_table(file, "/dl1/event/telescope/parameters/"+tel)['true_alt_tel'][0]

        except tables.exceptions.NoSuchNodeError:
            print(file, 'No pointing info!')
            continue

        time = read_table(file, "/dl1/event/telescope/parameters/"+tel)['local_time'][0]

        # find log 
        date = file.split('/')[-1].split('.')[0].split('_')[1]
        run = file.split('/')[-1].split('.')[0].split('_')[2]
        res = [i for i in logs if date in i]
        log_file = [i for i in res if run in i][0]

        # find fraction of pedestals in log
        word = 'Fraction of pedestal'
        with open(log_file, 'r') as fp:
            # read all lines in a list
            lines = fp.readlines()
            for line in lines:
                # check if string present on a current line
                if line.find(word) != -1:
                    ped_fraction = float(line.split(': ')[-1])
                    zenith_angles.append(zenith_angle)
                    times.append(time)
                    runs.append(int(run))
                    ped_fractions.append(ped_fraction)
                    filenames.append(file)

                    if "W1" in file:
                        wobble.append("W1")
                    elif "W2" in file:
                        wobble.append("W2")

    zenith_angles = np.array(zenith_angles)
    times = np.array(times)
    ped_fractions = np.array(ped_fractions)
    wobble = np.array(wobble)
    runs = np.array(runs)
    
    return zenith_angles, times, ped_fractions, wobble, runs, filenames


def moon_phase_angle(time):

    sun = get_sun(time)
    moon = get_body("moon", time)
    elongation = sun.separation(moon)
    return np.arctan2(sun.distance*np.sin(elongation),
                      moon.distance - sun.distance*np.cos(elongation))

def get_moon_params(data, config=None, tel=None, thinning=100):
    
    location = get_location(config=config, tel=tel)

    time = Time(data.iloc[::thinning, :]['local_time'], format='unix', scale='utc')
    moon = get_moon(time, location=location)
    horizon_frame = AltAz(obstime=time, location=location)
    moon_altaz = moon.transform_to(horizon_frame)
    
    moon_separation = ctaplot.ana.angular_separation_altaz(
    np.array(data.iloc[::thinning, :]['true_alt_tel']) * np.pi/180. * u.rad, 
    np.array(data.iloc[::thinning, :]['true_az_tel']) * np.pi/180. * u.rad, 
    moon_altaz.alt, 
    moon_altaz.az
    )
    
    return time, moon_altaz, moon_separation, moon_phase_angle(time)

def get_target_pos(target_name=None, ra=None, dec=None):
    """
    Get target coordinates as astropy SkyCoord. If target name
    is provided it tries to query the coordiantes. If it fails
    and manual coordinate input is not provided, then it returns
    SkyCoord with ra=0, dec=0.

    Parameters
    ----------
    target_name: string
    ra: float
        RA in degrees
    dec: float
        DEC in degrees
    Returns
    -------
    target_pos: astropy.coordinates.SkyCoord
    target_name: string

    """

    if target_name is None:
        logging.warning('Target name not provided! This is probably not a big issue..')
        target_name = 'UNKNOWN'
    try:
        target_pos = SkyCoord.from_name(target_name)
    except:
        logging.warning('Target not recognized. Manual coordinate input required.')

        if (ra is not None) and (dec is not None):
            target_pos = SkyCoord(ra  = ra *u.deg, dec = dec*u.deg, frame='icrs')
        else:
            logging.warning('Target coordinates not provided! Target coords set to (0,0) deg. This is probably not a big issue..')
            target_pos = SkyCoord(ra  = 0. *u.deg, dec = 0. *u.deg, frame='icrs')

    return target_pos, target_name


def get_GTIs(times):
    GTI_start = []
    GTI_stop = []

    t1 = times[0]

    GTI_start.append(t1)

    for time in np.array(times[1:-1]):
        t2 = time
        dt = t2-t1
        if dt>5:
            GTI_stop.append(t1)
            GTI_start.append(t2)
        t1 = time
    GTI_stop.append(times[-1])
    GTI = np.array([GTI_start*u.s,
                    GTI_stop *u.s])

    return GTI


def clip_alt(alt):
    """
    Make sure altitude is not larger than 90 deg 
    (it happens in some MC files for zenith=0),
    to keep astropy happy.

    This function is from lstchain
    https://github.com/cta-observatory/cta-lstchain
    Copyright (c) 2024 cta-lstchain developers
    Licensed under the 3-clause BSD style license.

    Parameters
    ----------
    alt: astropy.units.quantity.Quantity

    Returns
    -------
    astropy.units.quantity.Quantity

    """

    return np.clip(alt, -90.0 * u.deg, 90.0 * u.deg)


def add_disp_to_parameters(data):
    """
    We need to add DISP parameters to the dl1 table, 
    these are not implemented in ctapipe! This is a
    modification of the same functions from lstchain
    https://github.com/cta-observatory/cta-lstchain
    Copyright (c) 2024 cta-lstchain developers
    Licensed under the 3-clause BSD style license.

    Parameters
    ----------
    data: astropy.table.Table

    Returns
    -------
    tuple of astropy.units.quantity.Quantity

    """

    disp_parameters = disp(data['camera_frame_hillas_x'], # cog x
                            data['camera_frame_hillas_y'], # cog y
                            data['true_camera_x'],
                            data['true_camera_y'],
                            data['camera_frame_hillas_psi'].to(u.rad))

    return disp_parameters


def disp(
        cog_x, cog_y, src_x, 
        src_y, hillas_psi):
    """
    Calculates true DISP parameters shower 
    CoG and source coordinates. This is a
    modification of the same functions from lstchain
    https://github.com/cta-observatory/cta-lstchain
    Copyright (c) 2024 cta-lstchain developers
    Licensed under the 3-clause BSD style license.

    Parameters
    ----------
    cog_x: astropy.units.quantity.Quantity
        x coordinate of the shower center of gravity
    cog_y: astropy.units.quantity.Quantity
        y coordinate of the shower center of gravity
    src_x: astropy.units.quantity.Quantity
        x coordinate of the source position
        in the camera frame
    src_y: astropy.units.quantity.Quantity
        y coordinate of the source position
        in the camera frame
    hillas_psi: astropy.units.quantity.Quantity
        Angle between the x axis and the main axis
        of the shower

    Returns
    -------
    (disp_dx, disp_dy, disp_norm, disp_angle, disp_sign):
        disp_dx: astropy.units.quantity.Quantity
            Distance of the source from shower CoG
            along the x axis
        disp_dy: astropy.units.quantity.Quantity
            Distance of the source from shower CoG
            along the y axis
        disp_norm: astropy.units.quantity.Quantity
            Absolute value of a distance of the source 
            from shower CoG
        disp_angle: astropy.units.quantity.Quantity
            Hillas Psi angle
        disp_sign: numpy.ndarray
            Defines on which side of the main axis the 
            source lies

    """

    disp_dx = src_x - cog_x
    disp_dy = src_y - cog_y

    disp_norm = disp_dx * np.cos(hillas_psi) + disp_dy * np.sin(hillas_psi)
    disp_sign = np.sign(disp_norm)
    disp_norm = np.abs(disp_norm)

    # disp_sign : indicates in which direction, "positive" or "negative", we must move along the
    # reconstructed image axis (with direction defined by the versor cos(hillas_psi), sin(hillas_psi))
    # we must move from cog_x, cog_y to get closest to the true direction (src_x, src_y)

    if hasattr(disp_dx, '__len__'):
        disp_angle = np.arctan(disp_dy / disp_dx)
        disp_angle[disp_dx == 0] = np.pi / 2. * np.sign(disp_dy[disp_dx == 0])
    else:
        if disp_dx == 0:
            disp_angle = np.pi/2. * np.sign(disp_dy)
        else:
            disp_angle = np.arctan(disp_dy/disp_dx)

    return disp_dx, disp_dy, disp_norm, disp_angle, disp_sign


def disp_to_pos(disp_dx, disp_dy, cog_x, cog_y):
    """
    Calculates source coordinates from DISP params
    and CoG of the shower. This is a modification 
    of the same functions from lstchain
    https://github.com/cta-observatory/cta-lstchain
    Copyright (c) 2024 cta-lstchain developers
    Licensed under the 3-clause BSD style license.

    Parameters
    ----------
    disp_dx: astropy.units.quantity.Quantity
        Distance of the source from shower CoG
        along the x axis
    disp_dy: astropy.units.quantity.Quantity
        Distance of the source from shower CoG
        along the y axis
    cog_x: astropy.units.quantity.Quantity
        x coordinate of the shower center of gravity
    cog_y: astropy.units.quantity.Quantity
        y coordinate of the shower center of gravity

    Returns
    -------
    (source_pos_x, source_pos_y):
        source_pos_x: astropy.units.quantity.Quantity
        source_pos_y: astropy.units.quantity.Quantity

    """

    source_pos_x = cog_x + disp_dx
    source_pos_y = cog_y + disp_dy

    return source_pos_x, source_pos_y    


def disp_vector(disp_norm, disp_angle, disp_sign):
    """
    Calculates disp vector (source position relative
    to the position of the shower) from RF reconstructed
    quantities. This is a modification 
    of the same functions from lstchain
    https://github.com/cta-observatory/cta-lstchain
    Copyright (c) 2024 cta-lstchain developers
    Licensed under the 3-clause BSD style license.

    Parameters
    ----------
    disp_norm: numpy.ndarray
        Absolute value of a distance of the source 
        from shower CoG. Should be in metres
    disp_angle: pandas.core.series.Series
        Hillas Psi angle in rad
    disp_sign: numpy.ndarray
        Defines on which side of the main axis the 
        source lies

    Returns
    -------
    numpy.ndarray
        DISP vector with [disp_dx, disp_dy]

    """

    return np.transpose(polar_to_cartesian(disp_norm, disp_angle, disp_sign))


def polar_to_cartesian(norm, angle, sign):
    """
    Polar to cartesian transformation.
    As a convention, angle should be in 
    [-pi/2:pi/2]. This is a modification 
    of the same functions from lstchain
    https://github.com/cta-observatory/cta-lstchain
    Copyright (c) 2024 cta-lstchain developers
    Licensed under the 3-clause BSD style license.

    Parameters
    ----------
    norm: numpy.ndarray
    angle: numpy.ndarray
    sign: numpy.ndarray

    Returns
    -------
    (x, y):
        x: numpy.ndarray
        y: numpy.ndarray

    """

    assert np.isfinite([norm, angle, sign]).all()
    x = norm * sign * np.cos(angle)
    y = norm * sign * np.sin(angle)
    return x, y


def get_pointing_radec(input_file):
    """
    Gets array poiting in ra,dec from DL1 data file

    Parameters
    ----------
    input_file: string

    Returns
    -------
    pointing_radec: astropy.coordinates.SkyCoord

    """

    pointing = read_table(input_file, '/dl1/monitoring/subarray/pointing')
    p_ra = pointing['array_ra'][0] * pointing['array_ra'].unit
    p_dec = pointing['array_dec'][0] * pointing['array_dec'].unit
    pointing_radec = SkyCoord(ra=p_ra, dec=p_dec, frame='icrs')
    return pointing_radec


def get_subarray():

    subarray_file = pkg_resources.resource_filename(
                            'sst1mpipe',
                            path.join(
                                'data',
                                'sst1m_array.h5'
                            )
                        )
    subarray = SubarrayDescription.from_hdf(subarray_file, focal_length_choice="EQUIVALENT")
    return subarray


def get_cam_geom(tel=None):

    subarray = get_subarray()
    return subarray.tel[tel].camera.geometry


def get_dt_from_altaz(altaz_coord):

    subarray = get_subarray()

    d_tels = ((subarray.positions[21][1] - subarray.positions[22][1])**2 + 
              (subarray.positions[21][0] - subarray.positions[22][0])**2 +
              (subarray.positions[21][2] - subarray.positions[22][2])**2)**0.5

    tel_axis_az = (2*np.pi*u.rad-np.arctan(subarray.positions[22][1]/subarray.positions[22][0])).to('deg')
    tel_axis_alt = np.arcsin(2*subarray.positions[22][2]/d_tels)
    tel_axis_altaz = AltAz(alt=tel_axis_alt,az=tel_axis_az,
                           location=subarray.reference_location)
    DT = (d_tels/c.c * np.cos(tel_axis_altaz.separation( altaz_coord ))).to_value("ns")
    return DT
    

class Monitoring_R0_DL1:

    """
    Monitoring data for R0 -> DL1 processing

    Attributes
    ----------
        pointing_ra: float
            Pointing RA in deg
        pointing_dec: float
            Pointing DEC in deg
        target: string
            Observed source
        force_pointing: bool
            Pointing coordinates entered manualy 
            (i.e. not read from the raw FITS file)
        pointing_manual: bool
            Pointing coordinates entered manualy 
            (i.e. not read from the raw FITS file)
        wobble: string
        calibration_file: string
            Calibration file used
        window_file: string
            Window file used  
        n_saturated: int
            Total number of saturated events in the file
        n_pedestal: int
            Total number of pedestal events in the file
        n_survived_pedestals
            Number of pedestal events which survived cleaning
        n_triggered_tel1: int
            Total number of TEL1 triggered events in the file
        n_triggered_tel2: int
            Total number of TEL2 triggered events in the file
        swat_event_ids_used: bool
            If true, coincident events from given night have 
            exactly the same event ID
    """

    def __init__(self):
        self.pointing_ra = float
        self.pointing_dec = float
        self.force_pointing = bool
        self.pointing_manual = bool
        self.output_file = str
        self.output_file_px_charges = str
        self.window_file = str
        self.calibration_file = str
        self.target = str
        self.wobble = str
        self.n_pedestals = 0
        self.n_saturated = 0
        self.n_pedestals_survived = 0
        self.n_triggered_tel1 = 0
        self.n_triggered_tel2 = 0
        self.frac_rised = 0
        self.survived_charge_fraction_1 = []
        self.survived_charge_fraction_2 = []
        self.swat_event_ids_used = False


    def fill_target_info(self, input_file):

        target, ra_fits, dec_fits, wobble_fits = get_target(input_file, 
                                                            force_pointing=self.force_pointing
                                                            )
        self.target = target
        self.wobble = wobble_fits

        if (ra_fits is not None) & (dec_fits is not None):
            logging.info('Pointing info from the fits file: TARGET: ' + target + ', COORDS: ' + str(ra_fits) + ' ' + str(dec_fits) + ', WOBBLE: ' + wobble_fits)
            if self.force_pointing & (self.pointing_ra is not None) & (self.pointing_dec is not None):
                logging.info('Using pointing info from manual input anyway (forced by user).')
                self.pointing_manual=True
            else:
                self.pointing_ra = ra_fits
                self.pointing_dec = dec_fits
                self.output_file = self.output_file.split("_dl1.h5")[0] + "_" + wobble_fits + "_dl1.h5"
                self.output_file_px_charges = self.output_file_px_charges.split("_pedestal_hist.h5")[0] + "_" + wobble_fits + "_pedestal_hist.h5"
                self.pointing_manual=False

        elif force_pointing & (ra is not None) & (dec is not None):
            # Note that if there is no TARGET field in the fits file, it most probably means 
            # that it is a file where the shifters were tuning the trigger threshold on actual nsb conditions.
            logging.warning('No coordinates in the FITS header. Using pointing info from manual input. Are you sure that this is what you want?')
            self.pointing_manual=True

        else:
            logging.warning('No coordinates provided, exiting...')
            exit()


    def count_triggered(self, event, ismc=True):

        if ismc:
            tels = event.trigger.tels_with_trigger
            for tel in tels:
                if tel == 1:
                    self.n_triggered_tel1 += 1
                elif tel == 2:
                    self.n_triggered_tel2 += 1
        else:
            tel = event.sst1m.r0.tels_with_data[0]
            if tel == 21:
                self.n_triggered_tel1 += 1
            elif tel == 22:
                self.n_triggered_tel2 += 1

    def count_pedestals(self, event, ismc=True):

        if not ismc:
            tel = event.sst1m.r0.tels_with_data[0]
            event_type = event.sst1m.r0.tel[tel]._camera_event_type.value

            if event_type == 8:
                self.n_pedestals += 1

                if np.isfinite(event.dl1.tel[tel].parameters.hillas.intensity) :
                    self.n_pedestals_survived += 1

    def count_survived_charge(self, event, ismc=True):

        if ismc:
            tels = event.trigger.tels_with_trigger
            for tel in tels:
                cleaning_mask = event.dl1.tel[tel].image_mask
                if tel == 1:
                    self.survived_charge_fraction_1.append(sum(event.simulation.tel[tel].true_image[cleaning_mask])/sum(event.simulation.tel[tel].true_image))
                elif tel == 2:
                    self.survived_charge_fraction_2.append(sum(event.simulation.tel[tel].true_image[cleaning_mask])/sum(event.simulation.tel[tel].true_image))
                else:
                    logging.warning('Telescope %f not recognized, survived charge fraction not logged.', tel)


    def log_result_counts(self, ismc=True):

        logging.info('Total number of TEL1 triggered events in the file: %d', 
                    self.n_triggered_tel1)
        logging.info('Total number of TEL2 triggered events in the file: %d', 
                    self.n_triggered_tel2)
        
        if not ismc:
            logging.info('Total number of saturated events in the file: %d', 
                        self.n_saturated)
            logging.info('Total number of pedestal events in the file: %d', 
                        self.n_pedestals)
            if self.n_pedestals > 0:
                logging.info('Fraction of pedestal events that survived cleaning: %f', 
                            self.n_pedestals_survived/self.n_pedestals)
            else:
                logging.info('No pedestal events found!')

