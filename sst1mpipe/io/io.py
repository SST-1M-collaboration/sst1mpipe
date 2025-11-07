from ctapipe.io import read_table
from astropy.table import Table, QTable
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.table import join
import tables
import json
import os
from os import path
from traitlets.config import Config
import logging
from astropy.io.misc.hdf5 import (
    write_table_hdf5,
    read_table_hdf5
    )
from datetime import datetime

from sst1mpipe.utils import (
    add_pointing_mc,
    add_features,
    add_timing_features,
    add_log_true_energy,
    add_disp,
    event_selection,
    add_miss,
    get_telescopes,
    add_true_impact,
    get_location,
    add_event_id,
    get_tel_string,
    get_finite,
    get_pointing_radec,
    stereo_var_cuts
)

from astropy.io import fits

import sst1mpipe 

from sst1mpipe.io.containers import (
    DL1_info,
    DL2_info
)
from pyirf.cuts import (
    evaluate_binned_cut
)

import operator

from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
import glob
from gammapy.data import DataStore
import pkg_resources


def read_charges_data(file):
    """
    Read pixel charges distribution from 
    histogram files stored with sst1mpipe_r0_dl1
    from real data files

    Parameters
    ----------
    file: string
        Path to HDF file with histogram of
        pedestal charges

    Returns
    -------
    bin_centers: astropy.table.table.QTable
    charge_norm: astropy.table.table.Table

    """

    bins = read_table_hdf5(file, path='bins')['pe_bins']
    charge = read_table_hdf5(file, path='charge')['pixel_charge']
    charge_norm = charge/sum(charge)
    bin_centers = (bins[1:]-bins[:-1]) / 2 + bins[:-1]
    return bin_centers, charge_norm


def read_charges_mc(file, tel=None):
    """
    Read pixel charges distribution from 
    histogram files stored with sst1mpipe_r0_dl1
    from MC files (where we have distributions for 
    all simulated telescopes in single histogram
    file)

    Parameters
    ----------
    file: string
        Path to HDF file with histogram of
        pedestal charges
    tel: string

    Returns
    -------
    bin_centers: astropy.table.table.QTable
    charge_norm: astropy.table.table.Table

    """

    bins = read_table_hdf5(file, path='bins')['pe_bins']
    charge = read_table_hdf5(file, path='charge')['pixel_charge_tel'+str(tel)]
    charge_norm = charge/sum(charge)
    bin_centers = (bins[1:]-bins[:-1]) / 2 + bins[:-1]
    return bin_centers, charge_norm


def read_charge_images(file):

    ped_q_map = read_table_hdf5(file, path='charge_images')
    ped_q_map['start_time'] = np.array(ped_q_map['start_time'], dtype='datetime64[ns]')
    return ped_q_map


def load_config(cfg_file, ismc=False):
    """
    Reads sst1mpipe config file which must
    be stored as json file.

    Parameters
    ----------
    cfg_file: string
        Path to the config file
    ismc: bool

    Returns
    -------
    config: dict

    """
    if cfg_file is None:
        logging.info('No config file specified, loading default config.')
        if ismc:
            default_config = 'sst1mpipe_mc_config.json'
        else:
            default_config = 'sst1mpipe_data_config.json'
        cfg_file = pkg_resources.resource_filename('sst1mpipe',path.join('data',default_config))

    with open(cfg_file) as json_file:
            config = Config(json.load(json_file))

    return config


def load_source_catalog(catalog_file):
    """
    Reads sst1mpipe source catalog file which must
    be stored as json file.

    Parameters
    ----------
    catalog_file: string
        Path to the catalog file

    Returns
    -------
    catalog: dict

    """

    with open(catalog_file) as json_file:
            catalog = Config(json.load(json_file))
    return catalog


def check_outdir(outdir):
    """
    Checks whether the outpath exists and creates
    it if that is not the case.

    Parameters
    ----------
    outdir: string
        Path

    Returns
    -------

    """

    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
            logging.info("Created the output directory which didnt exist.")
        except FileExistsError:
            logging.info("Ah, the directory appeared in the meantime... \
                May happpen if you run more scripts in parallel terminals. \
                Should not be an issue.")


def write_dl2(
        dl2_table, output_file=None, telescope=None, 
        config=None, mode='w'):
    """
    Helper function to write DL2 file

    Parameters
    ----------
    dl2_table: pandas.DataFrame
    output_file: string
    telescope: string
    config: dict
    mode: string

    Returns
    -------

    """

    write_dl2_table(
        dl2_table, output_file=output_file, 
        table_path='/dl2/event/telescope/parameters', 
        table_name=telescope, config=config, 
        mode=mode
        )


def write_photon_list(
        pl_table, output_file=None, telescope=None, 
        config=None):
    """
    Helper function to write PL file

    Parameters
    ----------
    pl_table: pandas.DataFrame
        Photon list after gammaness cut
    output_file: string
    telescope: string
    config: dict

    Returns
    -------

    """

    write_pl_table(
        pl_table, output_file=output_file, 
        table_path='/photon_list', table_name=telescope, 
        config=config, mode='w'
        )


def write_pixel_charges_table(
        data, bin_edges, names=None, 
        output_file=None):

    res_t = Table(data=data, names=names)
    bins_t = QTable(data=bin_edges[..., np.newaxis], names=['pe_bins'])
    write_table_hdf5(bins_t, output_file, path='bins', overwrite=True, append=True, serialize_meta=True)
    write_table_hdf5(res_t, output_file, path='charge', overwrite=True, append=True)


def write_charge_images(ped_q_map, output_file=None):

    start_times = np.array(ped_q_map, dtype=object)[:,0]
    start_times = np.array([x.datetime64.astype(np.uint64) for x in start_times])[..., np.newaxis]
    start_times = np.hstack(start_times)

    n_ped = np.array(ped_q_map, dtype=object)[:,1]
    n_ped = np.hstack(n_ped.astype(int))

    q_map = np.array(ped_q_map, dtype=object)[:,2:3]
    q_map = np.hstack(q_map)
    q_map = np.vstack(q_map)

    stdev_map = np.array(ped_q_map, dtype=object)[:,3:4]
    stdev_map = np.hstack(stdev_map)
    stdev_map = np.vstack(stdev_map)

    ped_names = ['start_time', 'n', 'average_q', 'stdev_q']
    q_maps_t = Table(data=[start_times, n_ped, q_map, stdev_map], names=ped_names)
    write_table_hdf5(q_maps_t, output_file, path='charge_images', overwrite=True, append=True, serialize_meta=True)


def write_charge_fraction(file, survived_charge={}):
    """
    Writes fraction of survived charge after cleaning in
    the DL1 table of each telescope. Works on MC only
    as it needs true charge.

    Parameters
    ----------
    file: string
        Path to the output DL1 file
    survived_charge: dict
        Dictionary of survived charge fraction per
        telescope

    Returns
    -------

    """

    telescopes = get_telescopes(file)

    for tel in telescopes:

        params = read_table(file, "/dl1/event/telescope/parameters/" + tel)
        merged_table = params.copy()
        merged_table['survived_charge_fraction'] = survived_charge[tel]
        merged_table.write(file, path='/dl1/event/telescope/parameters/'+tel, overwrite=True, append=True) #, serialize_meta=True)


def write_extra_parameters(
        file, config=None, ismc=True, meanQ=None, wr_timestamps=None):
    """
    Opens the output DL1 file and adds some extra parameters
    to the DL1 table.

    Parameters
    ----------
    file: string
        Path to the output DL1 file
    config: dict
    ismc: bool
    meanQ: numpy.ndarray
        Mean charge from pedestal events
    wr_timestamps: numpy.array
        White Rabbit timestamps, where first colums is fulle seconds, and second is fraction seconds

    Returns
    -------

    """

    telescopes = get_telescopes(file)

    for tel in telescopes:

        params = read_table(file, "/dl1/event/telescope/parameters/" + tel)
        params = add_features(params)

        images = read_table(file, "/dl1/event/telescope/images/" + tel)
        params = add_timing_features(params, images)

        if ismc:
            params = add_pointing_mc(params, input_file=file, tel=tel)
            simulation = read_table(file, "/simulation/event/subarray/shower")
            merged_table = join(params, simulation, keys=['obs_id', 'event_id'])
            configuration = read_table(file, "/configuration/instrument/telescope/optics")
            merged_table['equivalent_focal_length'] = configuration['equivalent_focal_length'] * np.ones(len(merged_table))
            merged_table = add_disp(merged_table, config=config, telescope=tel)
            merged_table = add_miss(merged_table)
            merged_table = add_true_impact(merged_table, input_file=file, tel=tel)
            merged_table = add_log_true_energy(merged_table)

        else:
            
            # write WR timestamps
            if wr_timestamps is not None:
                params['time_wr_full_seconds'] = wr_timestamps[:, 0]
                params['time_wr_frac_seconds'] = wr_timestamps[:, 1]

            params["equivalent_focal_length"] = float(config['telescope_equivalent_focal_length'][tel])

            # Adding date
            merged_table = params.copy()
            date = int(file.split('/')[-1].split('_')[1])
            merged_table['date'] = date

            # adding trigger time
            trigger = read_table(file, path='dl1/event/subarray/trigger')

            # this is in seconds with numerical precision not enough for real stereo
            # Interestingly, this is most likely in UTC, i.e. local_time value stored in DL1 files
            # is by 37 sec smaller than white rabbit time
            merged_table['local_time'] = trigger['time'].unix

            #write mean charge from pedestal events to be used as NSB proxy for RF decision in DL1->DL2
            if meanQ is not None:
                meanQ = Table(data=meanQ, names=['obs_id', 'event_id', 'meanQ'], dtype=('i8', 'i8', 'f8'))
                merged_table = join(merged_table, meanQ, keys=['obs_id', 'event_id'])

        # Only combination of both overwrite and append works like expected, i.e. overwrite only telescope parameters and the rest 
        # of the DL1 file remains the same
        # NOTE: Unfortunately, we cannot store units with serialize_meta, because these are stored somehow weirdly as a new table and ctapipe merge tool
        # then cannot merge the files and raise error...
        merged_table.write(file, path='/dl1/event/telescope/parameters/'+tel, overwrite=True, append=True) #, serialize_meta=True)



def add_wr_dl1_stereo(file, dl1_data_tabs=[]):
    """
    Opens the DL1 stereo file after conincident event matching and
    adds two columns with high precision WR timestamps, which are
    stored in mono DL1 tables of both telescopes. This is neccessary
    because we use ctapipe DataWriter to store DL1 stereo file, but
    ctapipe containers ignore WR.

    Parameters
    ----------
    file: string
        Path
    dl1_data_tabs: list of pandas.DataFrame
        Mono DL1 tables of tel1 and tel2

    Returns
    -------

    """

    logging.info('Adding WR timestamps back into the DL1 stereo file...')
    telescopes = get_telescopes(file)

    if len(telescopes) == len(dl1_data_tabs):

        for tel in telescopes:

            params = read_table(file, "/dl1/event/telescope/parameters/" + tel)
            merged = params.copy()
            merged['time_wr_full_seconds'] = np.zeros(len(merged)).astype(np.int64)
            merged['time_wr_frac_seconds'] = np.zeros(len(merged)).astype(np.float64)
            merged['true_az_tel'] = np.zeros(len(merged)).astype(np.float64)
            merged['true_alt_tel'] = np.zeros(len(merged)).astype(np.float64)

            if tel == 'tel_021': params_tel = dl1_data_tabs[0]
            elif tel == 'tel_022': params_tel = dl1_data_tabs[1]

            # we cannot merge based on obs_id/event_id, because tel1/tel2 data has the same ids in the output file, but not in the input ones!
            # We also cannot merge based on only one parameter, because it turned out that the probability od having e.g. two events with the
            # very same (float64) intensity in data from a signle night is quite high
            params_tel_small = params_tel[['camera_frame_hillas_intensity', 'camera_frame_hillas_r', 'camera_frame_hillas_skewness', 'time_wr_full_seconds', 'time_wr_frac_seconds', 'true_az_tel', 'true_alt_tel']]

            for i, (intensity, r, skew) in enumerate(zip(np.array(params['camera_frame_hillas_intensity']), np.array(params['camera_frame_hillas_r']), np.array(params['camera_frame_hillas_skewness']))):
                mask = (intensity == params_tel_small['camera_frame_hillas_intensity']) & (r == params_tel_small['camera_frame_hillas_r']) & (skew == params_tel_small['camera_frame_hillas_skewness'])

                if sum(mask) == 1:
                    merged[i]['time_wr_full_seconds'] = params_tel_small[mask]['time_wr_full_seconds']
                    merged[i]['time_wr_frac_seconds'] = params_tel_small[mask]['time_wr_frac_seconds']
                    merged[i]['true_az_tel'] = params_tel_small[mask]['true_az_tel']
                    merged[i]['true_alt_tel'] = params_tel_small[mask]['true_alt_tel']

            merged.write(file, path='/dl1/event/telescope/parameters/'+tel, overwrite=True, append=True) #, serialize_meta=True)
            logging.info('WR timestamps added to ' + tel + ' param table.')

    else:
        logging.error('Different number of telescopes in the output file than number of provided DL1 tabs.')


def write_wr_timestamps(file, event_source=None):
    """
    Writes WR timestamps with high numerical precision as two 
    additional columns in the output DL1 table: time_wr_full_seconds,
    time_wr_frac_seconds. This is neccessary because the timestamp with 
    sufficient numerical precision can be extracted from event source only. 
    It is stored automaticaly in  dl1/event/subarray/trigger, from where 
    it can be read by ctapipe.io.read_table, but with low precision. 
    Therefore, if we want to have the timestamp in dl1, we need to read 
    it again from event source and store it in the existing dl1 file

    Parameters
    ----------
    file: string
        Path
    event_source: 
        sst1mpipe.io.sst1m_event_source.SST1MEventSource

    Returns
    -------

    """

    logging.info('Adding WR timestamps in DL1 table..')

    for i, event in enumerate(event_source):

        if i == 0:
            tel = event.sst1m.r0.tels_with_data[0]
            tel_string = get_tel_string(tel, mc=False)
            params = read_table(file, "/dl1/event/telescope/parameters/" + tel_string)
            time_wr_full_seconds_all = np.zeros(len(params)).astype(np.int64)
            time_wr_fractional_seconds_all = np.zeros(len(params)).astype(np.float64)

        event = add_event_id(event, filename=file, event_number=i)

        ev_mask = params['event_id'] == event.index.event_id

        if sum(ev_mask) == 1:

            localtime = event.sst1m.r0.tel[tel].local_camera_clock.astype(np.uint64)

            S_TO_NS = np.uint64(1e9)
            full_seconds = localtime // S_TO_NS
            fractional_seconds = (localtime % S_TO_NS) / S_TO_NS

            time_wr_full_seconds_all[ev_mask] = full_seconds
            time_wr_fractional_seconds_all[ev_mask] = fractional_seconds
    
    params['time_wr_full_seconds'] = time_wr_full_seconds_all
    params['time_wr_frac_seconds'] = time_wr_fractional_seconds_all

    params.write(file, path='/dl1/event/telescope/parameters/'+tel_string, overwrite=True, append=True) #, serialize_meta=True)
    logging.info('WR timestamps added to ' + tel_string + ' param table.')


def write_assumed_pointing(
        processing_info, config=None):
    """
    Writes pointing info (per event true_tel_az, true_tel_alt) 
    in the main DL1 table.

    Parameters
    ----------
    processing_info: 
        Class Monitoring_R0_DL1
    config: dict

    Returns
    -------

    """
    dl1_file = processing_info.output_file
    telescopes = get_telescopes(dl1_file)

    pointing_ra = float(processing_info.pointing_ra) * u.deg
    pointing_dec = float(processing_info.pointing_dec) * u.deg
    wobble_coords = SkyCoord(ra=pointing_ra, dec=pointing_dec, frame='icrs')

    for tel in telescopes:

        params = read_table(dl1_file, "/dl1/event/telescope/parameters/" + tel)
        location = get_location(config=config, tel=tel)
        time = Time(params['local_time'], format='unix', scale='utc')
        horizon_frame = AltAz(obstime=time, location=location)
        try:
            tel_pointing = wobble_coords.transform_to(horizon_frame)
            params['true_az_tel'] = tel_pointing.az.to_value(u.deg)
            params['true_alt_tel'] = tel_pointing.alt.to_value(u.deg)
        except ValueError:
            logging.info("Broken file", dl1_file)
            params['true_az_tel'] = np.nan
            params['true_alt_tel'] = np.nan

        # Only combination of both overwrite and append works like expected, i.e. overwrite only telescope parameters and the rest 
        # of the DL1 file remains the same
        # NOTE: Unfortunately, we cannot store units with serialize_meta, because these are stored somehow weirdly as a new table and ctapipe merge tool
        # then cannot merge the files and raise error...
        params.write(dl1_file, path='/dl1/event/telescope/parameters/'+tel, overwrite=True, append=True) #, serialize_meta=True)


def write_r1_dl1_cfg(file, config=None):
    """
    Write configuration of R1-D1 calibration
    in the output DL1 file.

    Parameters
    ----------
    file: string
        Path
    config: dict

    Returns
    -------

    """

    with tables.open_file(file, mode='a') as f:

        calibrator = config['CameraCalibrator']['image_extractor_type']
        image_processor = config['ImageProcessor']['image_cleaner_type']
        telescope_coords = config['telescope_coords']

        emin_cut_g = config['analysis']['gamma_min_simulated_energy_tev']
        emin_cut_p = config['analysis']['proton_min_simulated_energy_tev']

        for tel in telescope_coords.keys():

            t = f.create_table(
                '/configuration/r1_dl1/telescope_coords',
                tel,
                pd.DataFrame(config['telescope_coords'][tel], index=[0]).to_records(index=False),
                createparents=True,
            )

        t = f.create_table(
            '/configuration/r1_dl1/CameraCalibrator',
            calibrator,
            pd.DataFrame(config['CameraCalibrator'][calibrator], index=[0]).to_records(index=False),
            createparents=True,
        )

        t = f.create_table(
            '/configuration/r1_dl1/ImageProcessor',
            image_processor,
            pd.DataFrame(config['ImageProcessor'][image_processor], index=[0]).to_records(index=False),
            createparents=True,
        )

        t = f.create_table(
            '/configuration/r1_dl1',
            'Emin_cuts',
            pd.DataFrame({'emin_cut_g_tev': emin_cut_g, 'emin_cut_p_tev': emin_cut_p}, index=[0]).to_records(index=False),
            createparents=True,
        )


def load_more_dl1_tables_mono(
        file_list, config=None, check_finite=False, 
        time_min=0, time_max=np.inf, quality_cuts=False, 
        pointing_sel=None):
    """
    Loads and merges DL1 tables from multiple HDF DL1 mono files. 
    Tables are joined according to different criteria, such as
    time window or telescope pointing direction. 

    Parameters
    ----------
    file_list: list of strings
    config: dict
    check_finite: bool
        If True finite check of RF features is
        performed.
    time_min: float
        Adds DL1 table from given file to the merged table
        only if it at least partially overlaps with the 
        time interval specified by \'time_min\' and
        \'time_max\'. Unix time in seconds.
    time_max: float
        See \'time_min\'
    quality_cuts: bool
        If True event selection from the config 
        file is applied.
    pointing_sel: astropy.coordinates.SkyCoord
        Adds DL1 table from given file to the merged table
        only if the telescope pointing was within 0.2 deg 
        from the \'pointing_sel\' direction
     
    Returns
    -------
    dl1_data: pandas.DataFrame

    """

    good_files = 0
    bad_files = 0
    tel_file_list = []

    for dl1_file in file_list:

        tel = get_telescopes(dl1_file)
        
        if len(tel) == 1:

            df = load_dl1_sst1m(dl1_file, tel=tel[0], config=config, table='pandas', check_finite=check_finite, stereo=False, quality_cuts=quality_cuts)
            if len(df) > 0:

                if (min(df['local_time']) <= time_max) & (max(df['local_time']) >= time_min):

                    logging.info('Overlapping file with requested time interval found.')

                    if pointing_sel is not None:
                        pointing = get_pointing_radec(dl1_file)
                        if pointing.separation(pointing_sel) > 0.2 * u.deg:
                            logging.info('Skipping file with diferent pointing than requested.')
                            continue

                    if good_files == 0:
                        dl1_data = df
                        good_files += 1
                        tel_file_list.append(dl1_file)
                    else:
                        try:
                            dl1_data = pd.concat([dl1_data, df])
                            good_files += 1
                            tel_file_list.append(dl1_file)
                        except:
                            logging.warning("Skipping broken file: " + dl1_file)
                            bad_files += 1
                else:
                    logging.info('No overlap with requested time interval file, file skipped.')
            else:
                logging.info('No events left in the file, skipped.')
                bad_files += 1

        elif len(tel) == 0:
            logging.warning('No data in the file: ' + dl1_file)
            bad_files += 1
        else:
            logging.error('Data from more than one telescope already in the DL1 file! This is probably unexpected.')
            exit()

    logging.info('Data from ' + str(good_files) + ' good files merged.')
    logging.info('Bad files skipped:' + str(bad_files))
    logging.info('Total events: ' + str(len(dl1_data)))

    return dl1_data.reset_index()


def load_dl1_sst1m(
        input_file, tel=None, config=None, 
        table='astropy', check_finite=False, 
        stereo=False, quality_cuts=False, scale_intensities=False):
    """
    Reads DL1 table from the input HDF file.

    Parameters
    ----------
    input_file: string
        Path
    tel: string
        Either \'tel_00{1,2}\' (MC) or \'tel_02{1,2}\' (data)
    config: dict
    table: string
        Type of the output table: \'astropy\' 
        or \'pandas\'
    check_finite: bool
        If True finite check of RF features is
        performed.
    stereo: bool
        If True, extra columns with geometricaly reconstructed
        stereo quantities are added to the output table 
    quality_cuts: bool
        If True event selection from the config 
        file is applied.
    scale_intensities: bool
        If True a global scale on Hillas Intensity of each shower 
        is applied from config["NsbCalibrator"]["intensity_correction"]
    Returns
    -------
    data: pandas.DataFrame or astropy.table.Table

    """

    logging.info('Input file: %s', input_file) 
    events = read_table(input_file, "/dl1/event/telescope/parameters/" + tel)

    if scale_intensities:
        I_corr = float(config["NsbCalibrator"]["intensity_correction"][tel])
        logging.info('Intensity correction of %f applied on %s data.', I_corr, tel)
        events['camera_frame_hillas_intensity'] = events['camera_frame_hillas_intensity'] * I_corr

    if 'true_alt_tel' not in events.keys():
        try:
            pointing = read_table(input_file, "/dl1/monitoring/telescope/pointing/" + tel)
        except:
            logging.error('Adding pointing information failed! Pointing information is probably not stored in DL1 file.')
            exit()
        try:
            events['true_az_tel'] = pointing['azimuth'].to(u.deg).value
            events['true_alt_tel'] = pointing['altitude'].to(u.deg).value
        except:
            logging.error('Adding pointing information failed! Length of params and pointing tables probably dont match. Broken file.')
            exit()

    if stereo:
        stereo_impact = read_table(input_file, "/dl2/event/telescope/impact/HillasReconstructor/" + tel)
        stereo_geom = read_table(input_file, "/dl2/event/subarray/geometry/HillasReconstructor")
        merged = join(events.copy(), stereo_impact, keys=['obs_id', 'event_id', 'tel_id'], join_type='left', metadata_conflicts='silent')
        events = join(merged, stereo_geom, keys=['obs_id', 'event_id'], join_type='left', metadata_conflicts='silent')

    logging.info('Total N of events of %s: %d', tel, len(events))

    if quality_cuts & (config is not None):
        logging.info('Performing event selection.')
        events = event_selection(events, config=config)
    
    if check_finite & (config is not None):
        logging.info('Performing finite check.')
        events = get_finite(events, config=config, stereo=stereo)


    if stereo:
        logging.info('N of STEREO events of %s after selection cuts: %d', tel, len(events))
    else:
        logging.info('N of events of %s after selection cuts: %d', tel, len(events))

    if table == 'astropy':
    
        logging.info('DL1 data loaded as Astropy Table.')
        data = events

    elif table == 'pandas':

        logging.info('DL1 data loaded as Pandas Dataframe')

        # RFs like more pandas dataframe than astropy tables - it would be better to read the data in this format
        # at the first place. tels_with_trigger is multidimensional column and cannot be convereted to pandas df
        if 'tels_with_trigger' in events.keys():
            events.remove_column('tels_with_trigger')
        if stereo & ('HillasReconstructor_telescopes' in events.keys()):
            events.remove_column('HillasReconstructor_telescopes')
        data = events.to_pandas()

    return data


def load_dl2_sst1m(
        input_file, tel=None, config=None, 
        table='astropy', energy_min=0.0,
        scale_reco_energy=False):
    """
    Reads DL2 table from the input HDF file. If config
    is provided then event selection (according to the
    configuration) is applied. 

    Parameters
    ----------
    input_file: string
        Path
    tel: string
        Either \'tel_00{1,2}\' (MC) or \'tel_02{1,2}\' (data)
    config: dict
    table: string
        Type of the output table: \'astropy\' 
        or \'pandas\'
    energy_min: float
        Cut on minimum reconstructed energy [TeV]
    scale_reco_energy: bool
        If True a global scale on reconstructed energy of each shower 
        is applied from config["analysis"]["reco_energy_scaling_factor"]

    Returns
    -------
    data: pandas.DataFrame or astropy.table.Table

    """

    logging.info('Input file: %s', input_file) 
    events = read_table(input_file, "/dl2/event/telescope/parameters/" + tel)
    logging.info('Total N of events of %s: %d', tel, len(events))

    if config:
        logging.info('Performing event selection.')
        events = event_selection(events, config=config)
        logging.info('N of events of %s after selection cuts: %d', tel, len(events))

        if tel == 'stereo':
            events = stereo_var_cuts(events, config=config)


    if scale_reco_energy:
        logging.warning('Reconstructed energy scaled by %f', float(config["analysis"]["reco_energy_scaling_factor"]))
        events['reco_energy'] *= float(config["analysis"]["reco_energy_scaling_factor"])

    if energy_min > 0:
        events = events[events['reco_energy'] >= energy_min]
        logging.info('Cut on minumum reco energy > %f TeV applied.', energy_min)

    logging.info('N of DL2 events of %s after all input quality cuts: %d', tel, len(events))

    if table == 'astropy':
    
        logging.info('DL2 data loaded as Astropy Table.')
        data = events

    elif table == 'pandas':

        logging.info('DL2 data loaded as Pandas Dataframe')
        data = events.to_pandas()

    return data


def load_photon_list_sst1m(input_file, tel=None, config=None, table='astropy', energy_min=0.0):
    """
    Reads Photon list from the input HDF file. If config
    is provided then event selection (according to the
    configuration) is applied. 

    Parameters
    ----------
    input_file: string
        Path
    tel: string
        Either \'tel_00{1,2}\' (MC) or \'tel_02{1,2}\' (data)
    config: dict
    table: string
        Type of the output table: \'astropy\' 
        or \'pandas\'
    energy_min: float
        Cut on minimum reconstructed energy [TeV]

    Returns
    -------
    data: pandas.DataFrame or astropy.table.Table

    """

    logging.info('Input file: %s', input_file) 
    events = read_table(input_file, "/photon_list/" + tel)
    logging.info('Total N of events of %s: %d', tel, len(events))

    if config:
        logging.info('Performing event selection.')
        events = event_selection(events, config=config)

    if energy_min > 0:
        events = events[events['reco_energy'] >= energy_min]
        logging.info('Cut on minumum reco energy > %f TeV applied.', energy_min)

    logging.info('N of events of %s after selection cuts: %d', tel, len(events))

    if table == 'astropy':
    
        logging.info('Photon list loaded as Astropy Table.')
        data = events

    elif table == 'pandas':

        logging.info('Photon list loaded as Pandas Dataframe')
        data = events.to_pandas()

    return data


def load_dl1_pedestals(input_file):
    
    """
    Reads tables with pedestal info from the input HDF DL1 file.

    Parameters
    ----------
    input_file: string
        Path

    Returns
    -------
    pedestals: astropy.table.Table

    """

    pedestals = read_table(input_file, "/dl1/monitoring/telescope/pedestal")
    return pedestals


def write_dl1_pedestals(input_file, pedestal_table=None):

    """
    Write table of pedestal events from DL1 file into 
    another DL1 file. Typical usecase is to propagate 
    pedestals from mono DL1 to stereo DL1.

    Parameters
    ----------
    input_file: string
        Path
    
    pedestal_table: astropy.table.Table

    """

    try:
        write_table_hdf5(pedestal_table, input_file, 
            append=True, path='/dl1/monitoring/telescope/pedestal', 
            serialize_meta=False
            )
    except:
        logging.warning('Writing pedestals into the file failed!')


def load_extra_table(input_file, key=None, remove_column=None):

    table = read_table(input_file, key)
    if remove_column in table.keys():
        table.remove_column(remove_column)
    table_pd = table.to_pandas()

    return table_pd


def write_dl2_table(
        dl2, output_file=None, table_path=None, 
        table_name=None, config=None, mode='a'):
    """
    Opens the HDF file and writes DL2 table

    Parameters
    ----------
    dl2: pandas.DataFrame
    output_file: string
    table_path: string
    table_name: string
    config: dict
    mode: string

    Returns
    -------

    """

    for key in dl2.keys():

        if dl2[key].dtypes == object:
            logging.info('Removing key: %s, with content: %s', key, str(dl2[key].iloc[0]))
            logging.info('Because object saving is not supported in table')
            del dl2[key]
            continue
        
        # converting datatime columns in UNIX time in [ns]
        if dl2[key].dtypes == np.dtype('datetime64[ns]'):
            dl2[key] = dl2[key].astype('int')

    logging.info('Writing DL2 in %s', output_file)

    with tables.open_file(output_file, mode=mode) as f:

        full_node_path = f"{table_path}/{table_name}"
        if full_node_path in f:
            logging.info("Removing existing node %s", full_node_path)
            f.remove_node(table_path, name=table_name, recursive=True)

        logging.info("Creating new node %s", full_node_path)
        t = f.create_table(
            table_path,
            table_name,
            dl2.to_records(index=False),
            createparents=True,
        )

        if config:
            t.attrs["config"] = config


def write_pl_table(
        pl_table, output_file=None, table_path=None,
        table_name=None, config=None, mode='a'):
    """
    Opens the HDF file and writes PL table
    removing many keys from the DL2 table 
    which are not needed at this point.

    Parameters
    ----------
    pl_table: pandas.DataFrame
    output_file: string
    table_path: string
    table_name: string
    config: dict
    mode: string

    Returns
    -------

    """

    for key in pl_table.keys():

        if key in [
            'reco_disp_norm', 'reco_disp_sign', 'reco_disp_dx', 'reco_disp_dy',
            'camera_frame_hillas_skewness', 'camera_frame_hillas_kurtosis', 
            'camera_frame_hillas_x', 'camera_frame_hillas_y', 'camera_frame_hillas_r', 'camera_frame_hillas_phi',
            'camera_frame_hillas_length_uncertainty', 'camera_frame_hillas_width_uncertainty',
            'camera_frame_hillas_psi', 'camera_frame_timing_intercept', 'camera_frame_timing_deviation', 'camera_frame_timing_slope',
            'leakage_pixels_width_1', 'leakage_pixels_width_2', 'leakage_intensity_width_1',
            'concentration_cog', 'concentration_core', 'concentration_pixel', 'morphology_n_pixels',
            'morphology_n_small_islands', 'morphology_n_medium_islands', 'morphology_n_large_islands',
            'intensity_max', 'intensity_min', 'intensity_mean', 'intensity_std', 'intensity_skewness', 'intensity_kurtosis',
            'peak_time_max', 'peak_time_min', 'peak_time_mean', 'peak_time_std', 'peak_time_skewness', 'peak_time_kurtosis',
            'core_psi', 'log_camera_frame_hillas_intensity', 'camera_frame_hillas_wl'
        ]:
            logging.info('Removing key: %s.', key)
            del pl_table[key]
            continue

    logging.info('Writing Photon list in %s', output_file)

    with tables.open_file(output_file, mode=mode) as f:

        t = f.create_table(
            table_path,
            table_name,
            pl_table.to_records(index=False),
            createparents=True,
        )

        if config:
            t.attrs["config"] = config


def load_slow_data_bias_curve(file):
    """
    Reads slow control FITS file and extracts bias curve
    (rate scan) from it.

    Parameters
    ----------
    file: string

    Returns
    -------
    merge_data: dict

    """

    merge_data = {'timestamp': [], 'date': [],'biasCurveTriggerRate': [], 'biasCurveReadoutRate': [],
                  'biasCurvePatch7Threshold': [], 'biasCurveDroppedRate': [], 'appStatus': []}

    print("Loading file %s" % file)
    hdul = fits.open(file)
    bc_thr = hdul[1].data['biasCurvePatch7Threshold']
    bc_t = hdul[1].data['biasCurveTriggerRate']
    bc_r = hdul[1].data['biasCurveReadoutRate']
    bc_d = hdul[1].data['biasCurveDroppedRate']
    app_s = hdul[1].data['appStatus']
    timestamp = hdul[1].data['TIMESTAMP']

    thr_old = -1
    t_old = -1
    r_old = -1
    d_old = -1

    for i, vtime in enumerate(timestamp):
        if app_s[i][0] == 0 and bc_thr[i][0] > -0.5 and bc_t[i][0] > -0.5 and bc_d[i][0] > -0.5 and bc_r[i][0] > -0.5:
            if (t_old != bc_t[i][0] or thr_old != bc_thr[i][0] or r_old != bc_r[i][0] or d_old != bc_d[i][0]) and bc_t[i][0] > 0:
                thr_old = bc_thr[i][0]
                t_old = bc_t[i][0]
                r_old = bc_r[i][0]
                d_old = bc_d[i][0]
                merge_data['biasCurvePatch7Threshold'].append(bc_thr[i][bc_thr[i] > -0.5])
                merge_data['biasCurveTriggerRate'].append(bc_t[i][bc_t[i] > -0.5])
                merge_data['biasCurveReadoutRate'].append(bc_r[i][bc_r[i] > -0.5])
                merge_data['biasCurveDroppedRate'].append(bc_d[i][bc_d[i] > -0.5])
                merge_data['appStatus'].append(app_s[i][app_s[i] < 255])
                merge_data['timestamp'].append(vtime)
                merge_data['date'].append(datetime.fromtimestamp(float(vtime / 1000)))

    hdul.close()
    return merge_data


def load_drive_data(file):
    """
    Reads drive data FITS file and extracts
    timestamps from it.

    Parameters
    ----------
    file: string

    Returns
    -------
    merge_data: dict

    """

    merge_data = {'date': []}
    print('Opening file ' + file)
    # For the time being only reading splc states
    hdul = fits.open(file)
    cols = hdul[1].columns
    for j in cols.names:
        merge_data[j] = []
        merge_data[j].extend(hdul[1].data[j])
    merge_data['date'] = [datetime.fromtimestamp(float(k/1000)) for k in merge_data['TIMESTAMP']]
    hdul.close()
    return merge_data


def write_dl2_info(dl2_file, rfs_used=None):
    """
    Stores info tab with RF used and sst1mpipe 
    version in the DL2 file. These are necesary
    for DL3 production.

    Parameters
    ----------
    dl2_file: string
    rfs_used: string

    Returns
    -------

    """

    with tables.open_file(dl2_file, mode='a') as file:

        table = file.create_table(
            '/dl2',
            'info',
            DL2_info,
            "DL2 production info"
        )
        info = table.row
        info['sst1mpipe_version'] = sst1mpipe.__version__
        info['RF_used'] = rfs_used
        info.append()


def write_dl1_info(processing_info):
    """
    Stores info tab in the DL1 file

    Parameters
    ----------
        processing_info: Class Monitoring_R0_DL1

    Returns
    -------

    """
    dl1_file = processing_info.output_file

    with tables.open_file(dl1_file, mode='a') as file:

        table = file.create_table(
            '/dl1',
            'info',
            DL1_info,
            "DL1 production info"
        )
        info = table.row
        info['sst1mpipe_version'] = sst1mpipe.__version__
        info['target'] = processing_info.target
        info['ra'] = processing_info.pointing_ra
        info['dec'] = processing_info.pointing_dec
        info['manual_coords'] = processing_info.pointing_manual
        info['wobble'] = processing_info.wobble
        info['calib_file'] = processing_info.calibration_file
        info['window_file'] = processing_info.window_file
        info['n_saturated'] = processing_info.n_saturated
        info['n_pedestal'] = processing_info.n_pedestals
        info['n_survived_pedestals'] = processing_info.n_pedestals_survived
        info['n_triggered_tel1'] = processing_info.n_triggered_tel1
        info['n_triggered_tel2'] = processing_info.n_triggered_tel2
        info['swat_event_ids_used'] = processing_info.swat_event_ids_used
        info.append()


def get_dl1_info(file):

    info = read_table(file, "/dl1/info")
    return info


def load_more_dl2_files(files, config=None, 
                        gammaness_cut=None):
    """
    Loads and merges more DL2 hdf files. If config file is provided
    event selection is applied. Along with a single merged DL2 table
    the function returns also an array of Good Time Intervals, which
    can be used to define \'observation blocks\' (e.g. wobbles) for
    the purpose of DL3 creation.

    Parameters
    ----------
    files: list of strings
        List of paths to all dl2 files to be opened
    config: dict
    gammaness_cut: float or string
        Gammaness cut to be applied when loading and merging
        DL2 files (optional). If float, global gammaness cut
        is applied. Alternatively, it can be a string with
        path to the directory with gammaness cuts optimized 
        on MC. A set of gammaness cut for given zenith angle
        and NSB is found authomaticaly. The subdir structure 
        should follow the same logic as the RF model directories.

    Returns
    -------
    dl2_data: pandas.DataFrame
    times_all:
    GTI: numpy.ndarray
        Good Time Intervals. The splitting criteria are defined 
        as follows ():
            - pointing difference from previous DL2 file > 0.2 deg
            - time difference from previous event > 2 sec
            - RF used for for given DL2 file reconstruction is different
            from the previous DL2 file.
    """

    dl2_data = None
    pointing0 = None
    GTI_start, GTI_stop = [], []
    for input_file in files:

        try:
            tel_setup=get_telescopes(input_file, data_level="dl2")[0]
        except:
            logging.warning("No DL2 data in file : {}".format(input_file) )
            continue

        # Here we apply gammaness cut directly during the merge to save some memory.
        # local_times are merged without gammaness cut, which allows to split the observation
        # into different blocks later
        try:
            df0 = load_dl2_sst1m(input_file, tel=tel_setup, config=config, table='pandas', energy_min=0.0)

            # If there is no event in the file after selection cuts we skip it completely.
            if len(df0) == 0:
                logging.warning("No events in the file {} after selection cuts. SKIPPING (and not taking the time interval into account) ".format(input_file))
                continue

            times = df0['local_time']
            info = read_table(input_file,path='/dl2/info')
            RF_used = info['RF_used'][0].split('/')[-1]

            ## adding pointing info
            pointing = read_table(input_file,path='/dl1/monitoring/subarray/pointing')
            p_ra = np.unique(pointing['array_ra'])
            p_dec = np.unique(pointing['array_dec'])
            if (p_ra.shape[0]>1) or (p_dec.shape[0]>1):
                logging.warning('Multiple ra dec pointing in file {} -- we expect ony one!'.format(input_file))
            df0['array_ra']  = (p_ra[0]  * pointing['array_ra' ].unit).to_value('deg')
            df0['array_dec'] = (p_dec[0] * pointing['array_dec'].unit).to_value('deg')

            # GTI split based on pointing direction
            p_ra  = np.unique(df0['array_ra'])
            p_dec = np.unique(df0['array_dec'])
            pointing = SkyCoord(ra=p_ra * u.deg, dec=p_dec * u.deg, frame='icrs')

            # first file
            if pointing0 is None:
                pointing0 = pointing
                dt = 0
                GTI_start.append(times.iloc[0])
            else:
                dt = times.iloc[0] - times0.iloc[-1]

            # typical dt between events should be about 0.1 s, even for stereo it should not be higher that 0.5 s
            # length of a typical file is about 20 s, depending on conditions
            if (pointing.separation(pointing0) > 0.2 * u.deg) or (dt > 2):
                GTI_stop.append(times0.iloc[-1])
                GTI_start.append(times.iloc[0])
                pointing0 = pointing

            times0 = times

            if gammaness_cut is not None:
                if isfloat(gammaness_cut):
                    logging.info('Global gammaness cut {} applied.'.format(gammaness_cut))
                    mask = df0['gammaness'] > gammaness_cut
                    df = df0[mask].copy()
                    logging.info('N of events after gammaness cut: {}'.format(len(df)))
                else:
                    logging.info('Energy dependent gammaness cut applied.')

                    # Find gammaness file, there is only base path stored in the variable
                    if tel_setup[-1] == '1':
                        tel_mc = 'tel_001'
                    elif tel_setup[-1] == '2':
                        tel_mc = 'tel_002'
                    else:
                        tel_mc = tel_setup
                    try:
                        cut_file = glob.glob(gammaness_cut + '/' + RF_used + '/gammaness_cuts_*' + tel_mc + '*.h5')[0]
                        logging.info('Energy dendent cut table used: {}'.format(cut_file))
                        cut_table = read_table_hdf5(cut_file, path='gammaness_cuts')
                    except:
                        logging.warning("Cannot read gammaness cut file in the path: {}".format(cut_file))
                        cut_table = None

                    # This works only on pandas dataframe
                    df0_table = Table.from_pandas(df0)
                    mask_gg = evaluate_binned_cut(
                        df0_table["gammaness"],
                        df0_table["reco_energy"] * u.TeV,
                        cut_table,
                        operator.ge,
                    )
                    df = df0[mask_gg].copy()
                    logging.info('N of events after gammaness cut: {}'.format(len(df)))
            else:
                df = df0

            ## adding pipeline info
            df['RF_used'] = RF_used
            v1 = info['sst1mpipe_version'][0].split('.')[0]
            v2 = info['sst1mpipe_version'][0].split('.')[1]
            df['sst1mpipe_version'] = '{}_{}'.format(v1,v2)
            df['tel_setup'] = tel_setup
        except:
            logging.warning("Some problem with file {}".format(input_file))
            continue

        if dl2_data is None:
            dl2_data = df
            times_all = times
        else:
            try:
                dl2_data = pd.concat([dl2_data, df])
                times_all = pd.concat([times_all, times])
            except:
                logging.warning("Broken file", input_file)
                continue

    GTI_stop.append(times.iloc[-1])

    GTI = np.array([GTI_start * u.s,
                    GTI_stop * u.s])

    return dl2_data, times_all, GTI


def isfloat(num):
    try:
        float(num)
        return True
    except:
        return False


def load_distributions_sst1m(dist_path=None, dl3_path=None):
    """
    Reads all HDF files with DL1 distributions created with 
    sst1mpipe_extract_dl1_distributions per wobble and stack 
    them in output arrays.

    Parameters
    ----------
    dist_path: string
        Path to the HDF distribution files (one file per DL3)
    dl3_path: string
        Path to corresponding DL3 files.

    Returns
    -------
    histograms: numpy.ndarray
        Event rates bined in intensity
    histograms_diff: numpy.ndarray
        Differential event rates bined in intensity
    zeniths: numpy.ndarray
        Mean zenith angles per DL3 file
    obsids_sorted: numpy.ndarray
        OBS_ID of each DL3 file
    livetimes: numpy.ndarray
        Livetime of each DL3 file
    survived_ped: numpy.ndarray
        Fraction of survived pedestale events in each DL3 file
    bins: numpy.ndarray
        Bin edges of rate-intensity histograms
    """

    histograms = []
    histograms_diff = []
    zeniths = []
    obsids_sorted = []
    livetimes = []
    survived_ped = []

    no_ped = 0
    tables =  glob.glob(dist_path+'/intensity*.h5')
    if len(tables) == 0:
        # this is to read the data from date when we do not have pedestal events, so no recleaning, but we still index them in dl3 index files
        tables =  glob.glob(dist_path+'/intensity*.h5')
        no_ped = 1

    for table in tables:
        obsid = table.split('/')[-1].split('.')[0].split('_')[-1]
        hist = read_table(table, 'intensity_hist')
        t_elapsed = read_table(table, 't_elapsed')
        try:
            zenith = read_table(table, 'zenith')
        except:
            zenith = read_table(table, 'z_elapsed')

        # so that the pedestal fraction cut does not remove data for which we do not have pedestals (in those distribution files there is pedestal_frac=100.)
        if no_ped:
            survived_pedestal_frac = [0.]
        else:
            survived_pedestal_frac = read_table(table, 'survived_pedestal_frac')

        data_store = DataStore.from_dir(dl3_path)
        mask_datastore = data_store.obs_table['OBS_ID'] == int(obsid)
        if sum(mask_datastore) == 0:
            logging.warning('%d not in datastore!', int(obsid))
            continue

        histograms.append(np.array(hist['rate']).astype(np.float64))
        histograms_diff.append(np.array(hist['diff_rate']).astype(np.float64))
        zeniths.append(np.array(zenith)[0].astype(np.float64))
        obsids_sorted.append(int(obsid))
        survived_ped.append(np.array(survived_pedestal_frac)[0].astype(np.float64))
        bins = np.hstack((np.array(hist['low']),hist['high'][-1]))    
        mask_datastore = data_store.obs_table['OBS_ID'] == int(obsid)
        livetime = data_store.obs_table['LIVETIME'].value
        livetimes.append(livetime[mask_datastore][0])
    
    histograms = np.array(histograms)
    histograms_diff = np.array(histograms_diff)
    zeniths = np.array(zeniths)
    survived_ped = np.array(survived_ped)
    obsids_sorted = np.array(obsids_sorted)
    livetimes = np.array(livetimes).flatten()

    return histograms, histograms_diff, zeniths, obsids_sorted, livetimes, survived_ped, bins


def get_used_qe_simtel(source):
    """
    Reads the simtel cfg file used to run current MC production and extracts
    all Quntum efficiency (PDE) files listed there. There are also default (dummy)
    PDE files there as the proper fields in simtel must be first initialized with
    something before running. This is not dangerous as long us one does not 
    use one of the real PDE files used for the real production as the dummy file..

    Parameters
    ----------
    source: ctapipe.io.EventSource

    Returns
    -------
    used_qe: numpy.array
        Array of used PDEs to produce given simtel file

    """

    used_qe = []
    stringlist=[x[1].decode('utf-8') for x in source.file_.history]
    string_array = np.array(stringlist)
    indices = np.arange(0, len(string_array))
    index_array_eff = indices[np.char.find(string_array, 'QUANTUM_EFFiciency') != -1]
    for string in string_array[index_array_eff]:
        s = string.replace('QUANTUM_EFFiciency ', '').split('%')[0].strip()
        used_qe.append(s[:s.rfind('.')])
    return np.array(used_qe)


def get_pde_correction_factors():
    """
    Reads the default calibration file containing the PDE
    corrections for different PDE files used in MC production.
    The file is expected to be stored in ../data/mc_pde_correction_factors.json

    Returns
    -------
    pde_corr: dict

    """

    try:
        default_pde_corr_file = 'mc_pde_correction_factors.json'
        pde_corr_file = pkg_resources.resource_filename('sst1mpipe',path.join('data',default_pde_corr_file))

        with open(pde_corr_file) as json_file:
                pde_corr = Config(json.load(json_file))
        return pde_corr
    except:
        logging.error('%s file not found! Either make sure it is there, or turn off the PDE correction in the cfg file!'.format(pde_corr_file))
        exit()