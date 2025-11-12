from astropy.table import Table
import numpy as np
import pandas as pd
import logging
import os
import glob

from sst1mpipe.utils import (
    get_telescopes,
    get_horizon_frame,
    get_GTIs,
)
from sst1mpipe.analysis import add_reco_ra_dec
from sst1mpipe.io import load_config, load_more_dl2_files, check_outdir

from ctapipe.coordinates import TelescopeFrame

import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
from astropy.io.misc.hdf5 import read_table_hdf5

from gammapy.data import DataStore


def photon_df_to_fits(dl2_photons,
                      obs_id = None,
                      target_name = None,
                      target_pos = None,
                      config = None,
                      RF_used = None,
                      start_t = None,
                      end_t   = None,
                      gammaness_cuts = None,
                      event_class    = None,
                      ):
    """
    Create HDU list to be stored in the DL3 fits file
    from input photon list.

    Parameters
    ----------
    dl2_photons: pandas.DataFrame
        Merged DL2 table with gammaness cut applied 
        (photon list) for one single GTI
    obs_id: float
        Usualy obs_id of the first DL2 file in given 
        bunch after the merge
    target_name: string
    target_pos: astropy.coordinates.SkyCoord
    config: dict
    RF_used: string
    start_t: float
        GTI time bin edge
    end_t: float
        GTI time bin edge
    gammaness_cuts: None or string
        Flags the IRF assigned to 
        given DL3 file. If None, global gammaness 
        cut from the config file is applied.
    event_class: int
        Event class in given photon list to find a proper IRF
    Returns
    -------
    hdu_list: 
        list of astropy.io.fits.hdu.table.BinTableHDU

    """

    tel_setup = dl2_photons['tel_setup'].iloc[0]
    instrument = "SST1M_"+tel_setup
    if tel_setup=='stereo':
        NTEL=2
    else:
        NTEL =1

    if start_t is None :
        start_t = dl2_photons['local_time'].min()

    if end_t is None :
        end_t = dl2_photons['local_time'].max()

    GTI = np.array([[start_t]*u.s,[end_t]*u.s])
    
    ## Livetime
    livetime = end_t-start_t

    ### ra dec array pointing ###
    p_ra  = np.unique(dl2_photons['array_ra'])
    p_dec = np.unique(dl2_photons['array_dec'])
    if (p_ra.shape[0]>1) or (p_dec.shape[0]>1):
        logging.error('Multiple ra dec pointing in DF -- we expect only one!!!! (obs_id {})'.format(obs_id))

    radec_pointing = SkyCoord(ra  = p_ra[0] *u.deg,
                              dec = p_dec[0]*u.deg,
                              frame='icrs')

    ## reco direction in FOV
    tel_pointing = AltAz(alt = np.array(dl2_photons['true_alt_tel']) * u.deg,
                         az  = np.array(dl2_photons['true_az_tel'])*u.deg)
    mean_zen = 90 - dl2_photons['true_alt_tel'].mean()

    telescopeframe = TelescopeFrame(telescope_pointing=tel_pointing)
    reco_altaz = AltAz(alt= np.array(dl2_photons['reco_alt'])*u.deg,
                       az = np.array(dl2_photons['reco_az'])*u.deg)

    offset = tel_pointing.separation(AltAz(alt = np.array(dl2_photons['reco_alt'])*u.deg,
                                           az =  np.array(dl2_photons['reco_az'])*u.deg))

    if ('reco_ra' not in dl2_photons.keys()) or ('reco_dec' not in dl2_photons.keys()):
        horizon_frame = get_horizon_frame(config=config, telescope='tel_021', times=Time(dl2_photons['local_time'], format='unix'))
        dl2_photons = add_reco_ra_dec(dl2_photons, horizon_frame=horizon_frame)
        logging.info('Calculated reco RA and DEC, which were missing in the photon list')

    ## IRF name
    #zmean = 90 - dl2_photons['true_alt_tel'].mean()
    if gammaness_cuts is None:
        gammaness_cut = config['analysis']['global_gammaness_cut']
    else: 
        gammaness_cut = 'energydep'

    # Removing NSB info from the RF_used string
    try:
        name = RF_used.split('_nsb')[0]
    except:
        name = RF_used
    if event_class is not None:
        ec_str = "_ec{}".format(event_class)
    else: ec_str = ""
    irf_name = "{}_gc{}{}".format(name, gammaness_cut, ec_str)
    logging.info('Expected IRF name for given DL3 file: %s', irf_name)
    pipeline_version = dl2_photons["sst1mpipe_version"].iloc[0]

    ########### Build TABLE
    dl3_columns =["EVENT_ID","TIME","RA","DEC","ENERGY","DETX","DETY"]
    T_event = Table()
    T_event.add_column(dl2_photons['event_id'],name="EVENT_ID")
    T_event.add_column(np.array(dl2_photons['local_time'])*u.s,name="TIME")

    T_event.add_column(np.array(dl2_photons['reco_ra'])*u.deg,name="RA")
    T_event.add_column(np.array(dl2_photons['reco_dec'])*u.deg,name="DEC")

    T_event.add_column(np.array(dl2_photons['reco_energy'])*u.TeV,name="ENERGY")

    T_event.add_column(reco_altaz.transform_to(telescopeframe).fov_lat,name="DETX")
    T_event.add_column(reco_altaz.transform_to(telescopeframe).fov_lon,name="DETY")
    fbt_events = fits.BinTableHDU(T_event,name='EVENTS')

    fbt_events.header['CREATOR']  = ('sst1mpipe_data_dl2_dl3','Program which created the file')             
    fbt_events.header['TELESCOP'] = (instrument,'Telescope')                                      
    fbt_events.header['OBS_ID']   = (obs_id, 'Observation identifier')                         
    fbt_events.header['DATE_OBS'] = (Time(start_t,
                                          format = 'unix',
                                          scale  = 'tai').strftime('%Y-%m-%d'),
                                     'Observation start date')                         
    fbt_events.header['TIME_OBS'] = (Time(start_t,
                                          format = 'unix',
                                          scale  = 'tai').strftime('%H:%M:%S'),
                                     'Observation start time')

    fbt_events.header['DATE_END'] = (Time(end_t,
                                          format = 'unix',
                                          scale  = 'tai').strftime('%Y-%m-%d'),
                                     'Observation end date')                       
    fbt_events.header['TIME_END'] = (Time(end_t,
                                          format = 'unix',
                                          scale  = 'tai').strftime('%H:%M:%S'),
                                     'Observation end time')          

    fbt_events.header['TSTART']   = (start_t,'[s] Mission time of start of observation')
    fbt_events.header['TSTOP']    = (end_t,'[s] Mission time of end of observation')

    fbt_events.header['MJDREFI']  = (40587,'[days] Integer part of time reference MJD') ## UNIX TIME REF 
    fbt_events.header['MJDREFF']  = (0,'[days] Fractional part of time reference MJD' )  
    fbt_events.header['TIMEUNIT'] = ('s',     'Time unit')                                     
    fbt_events.header['TIMESYS']  = ('TAI',   'Time system')                                    
    fbt_events.header['TIMEREF']  = ('TOPOCENTER', 'Time reference')

    fbt_events.header['TELAPSE']  = (livetime, '[s] Mission elapsed time')                       
    fbt_events.header['ONTIME']   = (livetime, '[s] Total good time including deadtime')     
    fbt_events.header['LIVETIME'] = (livetime, '[s] Total livetime')                          
    fbt_events.header['DEADC']    = (1, 'Deadtime correction factor') ## ASSUMING NO DEADTIME
    fbt_events.header['TIMEDEL']  = (1, 'Time resolution')

    fbt_events.header['OBJECT']   = (target_name, 'Observed object')                             
    fbt_events.header['RA_OBJ']   = (target_pos.ra.to_value('deg'), '[deg] Target Right Ascension')                
    fbt_events.header['DEC_OBJ']  = (target_pos.dec.to_value('deg'), '[deg] Target Declination')

    fbt_events.header['RA_PNT']   = (radec_pointing.ra.to_value("deg"),  '[deg] Pointing Right Ascension')            
    fbt_events.header['DEC_PNT']  = (radec_pointing.dec.to_value("deg"), '[deg] Pointing Declination')

    fbt_events.header['ALT_PNT']  = (dl2_photons['true_alt_tel'].mean(), '[deg] Average altitude of pointing')       
    fbt_events.header['AZ_PNT']   = (dl2_photons['true_az_tel'].mean(),  '[deg] Average azimuth of pointing ')             
    fbt_events.header['RADECSYS'] = ('FK5','Coordinate system')                              
    fbt_events.header['EQUINOX']  = (2000,'Epoch')                                          
    fbt_events.header['CONV_DEP'] = (0, 'Convergence depth of telescopes')               
    fbt_events.header['CONV_RA']  = (0, '[deg] Convergence Right Ascension')              
    fbt_events.header['CONV_DEC'] = (0, '[deg] Convergence Declinason')   
    fbt_events.header['OBSERVER'] = ('SST1M Collaboration','Observer')                                       
    fbt_events.header['N_TELS']   = (NTEL, 'Number of telescopes in event list')             
    fbt_events.header['TELLIST']  = ('21-22','Telescope IDs')
    ## subarray.reference_location.to_geodetic()
    ## lon=<Longitude 14.782501 deg>, lat=<Latitude 49.912442 deg>, height=<Quantity 511.195 m>)
    fbt_events.header['GEOLAT']  =    (config['array_center_coords']['lat_deg'], '[deg] Geographic latitude of array centre')
    fbt_events.header['GEOLON']  =    (config['array_center_coords']['long_deg'], '[deg] Geographic longitude of array centre')    
    fbt_events.header['ALTITUDE']=    (config['array_center_coords']['height_m'],  '[km] Altitude of array centre')

    fbt_events.header['EUNIT']   = ('TeV',           'Energy unit')                                    
    fbt_events.header['EVTVER']  = (pipeline_version,        'Event list version number')                      
    fbt_events.header['CALDB']   = (pipeline_version,             'Calibration database')                           
    fbt_events.header['IRF']     = (irf_name, 'Instrument Response Function')

    fbt_gti = GTI_df_to_hdu(GTI)
    hdu_list = [fits.PrimaryHDU(), fbt_events, fbt_gti]

    return hdu_list


def GTI_df_to_hdu(GTI):
    """
    Creates HDU table from GTI array containing 
    time bin edges for each observation block

    Parameters
    ----------
    GTI: np.ndarray

    Returns
    -------
    fbt_gti: astropy.io.fits.hdu.table.BinTableHDU

    """

    T_GTI = Table()
    T_GTI.add_column(GTI[0]*u.s,name='START')
    T_GTI.add_column(GTI[1]*u.s,name='STOP')

    fbt_gti = fits.BinTableHDU(T_GTI,name="GTI")
    ## UNIX TIME T0
    fbt_gti.header['MJDREFI']  = (40587,'[days] Integer part of time reference MJD')      
    fbt_gti.header['MJDREFF']  = (0 ,'[days] Fractional part of time reference MJD' )  
    fbt_gti.header['TIMEUNIT'] = ('s',     'Time unit')                                     
    fbt_gti.header['TIMESYS']  = ('TAI',   'Time system')            ## TBC
    fbt_gti.header['TIMEREF']  = ('TOPOCENTER', 'Time reference')    ## TBC
    print(type(fbt_gti))
    return fbt_gti


def create_hdu_index(file_list, out_dir=None, irf_dir=None):
    """
    Create index files linking DL3 fits files with IRFs.

    Parameters
    ----------
    file_list: list
        DL3 file list
    out_dir: string
    irf_dir: string
        Top level directory with IRFs. Specific IRF file in 
        ../data/ is found authomaticaly to all given observing 
        blocks.

    Returns
    -------

    """

    os.environ['CALDB'] = irf_dir

    try:
        data_store = DataStore.from_events_files(file_list)
    except:
        logging.error("Likely, No IRFs exists for those DL3 file... We cant create HDU-index")
        exit()
    ## Hack to load proper name of the irf in hdu tables (for now hard coded in gammapy)
    for ii in np.where(data_store.hdu_table['HDU_TYPE']=='psf')[0]:
        data_store.hdu_table[ii]["HDU_CLASS"]= "psf_table"
    for ii in np.where(data_store.hdu_table['HDU_TYPE']=='bkg')[0]:
        data_store.hdu_table[ii]["HDU_CLASS"]= "bkg_2d"
 
    ## writing files 
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(data_store.hdu_table)
    ]).writeto(out_dir+"/hdu-index.fits.gz", overwrite=True)

    fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(data_store.obs_table)
    ]).writeto(out_dir+"/obs-index.fits.gz", overwrite=True)


def dl2_dir_to_dl3(target_name   = None,
                   target_pos    = None,
                   dl2_dir       = None,
                   config_file   = None,
                   out_dir       = None,
                   gammaness_cuts= None,
                   event_classes   = None,
                   ):
    """
    Process all DL2 files in given directory in DL3. DL2 files
    are merged into several bunches indexed by \'obs_id\'. 
    The clustering is driven by GTIs and mostly result in
    one DL3 file per wobble, but there are also other criteria
    such as RF used (wobble can be further split if the pointing
    passed through an edge of RF zenith angle bin).

    Parameters
    ----------
    target_name: string
    target_pos: astropy.coordinates.SkyCoord
    dl2_dir: str
        Path to target DL2 directory in which the files
        are to be processed to DL3.
    config_file: dict
    out_dir: string
    gammaness_cuts: string
        Gammaness cut to be applied when loading and merging
        DL2 files (optional). If None, global gammaness cut
        is applied. Alternatively, it can be a string with
        path to the directory with gammaness cuts optimized 
        on MC. A set of gammaness cut for given zenith angle
        and NSB is found authomaticaly. The subdir structure 
        should follow the same logic as the RF model directories.
    event_classes: list of ints
        classes in DL2 column event_type to be processed.
    Returns
    -------
    created_files: list of paths

    """
    
    config = load_config(config_file)

    created_files = {}
    created_files['all'] = []

    all_files = glob.glob(dl2_dir+'/*dl2.h5')
    all_files.sort()

    if len(all_files)==0:
        return []

    if gammaness_cuts is None:
        df_dl2_unsort, times_unsort, GTIs = load_more_dl2_files(
                all_files, 
                config = config, 
                gammaness_cut=config['analysis']['global_gammaness_cut']
                )
    else:
        df_dl2_unsort, times_unsort, GTIs = load_more_dl2_files(
                all_files, 
                config = config, 
                gammaness_cut=gammaness_cuts
                )
    df_dl2 = df_dl2_unsort.copy().sort_values("local_time").reset_index()
    #times = np.array(times_unsort)
    #times.sort()

    # check outdirs
    if len(event_classes):
        for evtclass in event_classes:
            outdir = out_dir + '/event_class_' + str(evtclass)
            check_outdir(outdir)

    logging.info("{} GTIs (i.e. wobbles) in total.".format(len(GTIs.T)-1))
    
    for obs_ii, GTI in enumerate(GTIs.T):

        df_t = df_dl2[(df_dl2['local_time']>GTI[0]) &
                      (df_dl2['local_time']<GTI[1]) ]
        if np.array(df_t).shape[0]==0:
            logging.warning("Skipping empty df {}".format(obs_ii)) ## MORE INFO ? 
            continue

        for RF_used in np.unique(df_t['RF_used']):
            df_tt = df_t[df_t['RF_used']==RF_used]

            ## TODO better handle GTI ? how? ... 
            if np.unique(df_t['RF_used']).shape[0]>1:
                start_t = df_tt['local_time'].min()
                end_t   = df_tt['local_time'].max()
            else:
                start_t = GTI[0]
                end_t   = GTI[1]

            # Not needed, this event selection is performed when dl2 data are loaded. 
            #df_tt = event_selection(df_tt, config=config)
            
            if len(df_tt)<1:
                continue

            ## better deefinition of obs_id ?
            obs_id = df_tt["obs_id"].min()

            if len(event_classes):
                for evtclass in event_classes:
                    out_dir_c = out_dir + '/event_class_' + str(evtclass)

                    logging.info("Making DL3s for event class {}".format(evtclass))
                    df_tt_c = df_tt[df_tt['event_type'] == evtclass]
                    logging.info("{} events of class {}".format(len(df_tt_c), evtclass))

                    # If there are no photons, the DL3 file is not created. This may potentialy
                    # raise some issues with exposure not being properly calculated.
                    if len(df_tt_c):
                        hdulist = photon_df_to_fits(df_tt_c,
                                                    obs_id           = obs_id,
                                                    target_name      = target_name,
                                                    target_pos       = target_pos,
                                                    config           = config,
                                                    RF_used          = RF_used,
                                                    start_t          = start_t,
                                                    end_t            = end_t,
                                                    gammaness_cuts   = gammaness_cuts,
                                                    event_class      = evtclass,
                                                    )
                        outname = os.path.join(out_dir_c,"SST1M_{}_obs_id_{}_dl3.fits".format(target_name,obs_id))
                        fits.HDUList(hdulist).writeto(outname, overwrite=True)
                        key = 'event_class_'+str(evtclass)
                        if key not in created_files:
                            created_files[key] = []
                        created_files[key].append(outname)
                    else:
                        outname = os.path.join(out_dir_c,"SST1M_{}_obs_id_{}_dl3.fits".format(target_name,obs_id))
                        logging.warning("No events of class {}, file {} not created!".format(evtclass,outname))
            else:
                if len(df_tt):
                    hdulist = photon_df_to_fits(df_tt,
                                                    obs_id           = obs_id,
                                                    target_name      = target_name,
                                                    target_pos       = target_pos,
                                                    config           = config,
                                                    RF_used          = RF_used,
                                                    start_t          = start_t,
                                                    end_t            = end_t,
                                                    gammaness_cuts = gammaness_cuts
                                                    )
                    outname = os.path.join(out_dir,"SST1M_{}_obs_id_{}_dl3.fits".format(target_name,obs_id))
                    fits.HDUList(hdulist).writeto(outname, overwrite=True)
                    created_files['all'].append(outname)
                else:
                    outname = os.path.join(out_dir,"SST1M_{}_obs_id_{}_dl3.fits".format(target_name,obs_id))
                    logging.warning("No events for this wobble/RF, file {} not created!".format(outname))

    return created_files
