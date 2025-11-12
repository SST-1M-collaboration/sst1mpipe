from astropy.coordinates import (
    SkyCoord, 
    EarthLocation, 
    AltAz, 
    Angle, 
    SkyOffsetFrame
)
from astropy.time import Time
from ctapipe.coordinates import CameraFrame
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from gammapy.stats import WStatCountsStatistic
import math

from sst1mpipe.utils import (
    clip_alt,
    get_horizon_frame,
)

from sst1mpipe.performance import get_theta
import logging
from gammapy.data import DataStore


def add_reco_ra_dec(data, horizon_frame=None):
    """
    Calculate reconstructed RA, DEC from reconstructed ALT,AZ in DL2 table 
    for timestamps specified in horizon_frame and add them as new columns

    Parameters
    ----------
    data: pandas.DataFrame with DL2 table
    horizon_frame: astropy.coordinates.AltAz

    Returns
    -------
    pandas.DataFrame

    """

    data = data.copy()
    reco_coord = SkyCoord(
        alt=np.array(data['reco_alt'])*u.deg, 
        az=np.array(data['reco_az'])*u.deg, 
        frame=horizon_frame
        )
    radec = reco_coord.transform_to('icrs')
    if 'reco_ra' in data.columns:
        logging.warning('Reconstructed RA DEC already in the data, replacing..')
        data = data.drop(columns=['reco_ra', 'reco_dec']).copy()
    data['reco_ra'] = radec.ra.value
    data['reco_dec'] = radec.dec.value

    return data


def add_source_altaz(data, source=None, horizon_frame=None):
    """
    Add altitude and azimuth of the source position to each event
    in input DL2 table

    Parameters
    ----------
    data: pandas.DataFrame with DL2 table
    source: astropy.coordinates.SkyCoord
    horizon_frame: astropy.coordinates.AltAz

    Returns
    -------
    pandas.DataFrame

    """

    data = data.copy()

    source_altaz = source.transform_to(horizon_frame)
    if 'source_alt' in data.columns:
        logging.warning('Source coords already in the data, replacing..')
        data = data.drop(columns=['source_alt', 'source_az']).copy()
    data['source_alt'] = source_altaz.alt.value
    data['source_az'] = source_altaz.az.value

    return data


def get_camera_frame(data, config=None, telescope=None):
    """

    Parameters
    ----------
    data: pandas.DataFrame with DL1 or DL2 table
    config: dict
    telescope: str (tel_021/tel_022)

    Returns
    -------
    ctapipe.coordinates.CameraFrame

    """

    horizon_frame = get_horizon_frame(
        config=config, 
        telescope=telescope, 
        times=Time(data.local_time, format='unix', scale='utc')
        )
    focal = np.array(data["equivalent_focal_length"])[0] * u.m
    logging.info('Focal length used: {}'.format(np.array(data["equivalent_focal_length"])[0]))
    ponting_coords = SkyCoord(
        alt = np.array(data['true_alt_tel'])*u.deg, 
        az = np.array(data['true_az_tel'])*u.deg, 
        frame=horizon_frame
        )

    camera_frame = CameraFrame(
        focal_length=focal, 
        telescope_pointing=ponting_coords, 
        obstime=horizon_frame.obstime,
        location=horizon_frame.location)

    return camera_frame


def add_source_xy(data, source=None, camera_frame=None):
    """
    Add source x,y coordinates in camera frame to each event
    in input dataframe

    Parameters
    ----------
    data: pandas.DataFrame with DL2 table
    source: astropy.coordinates.SkyCoord
    camera_frame: ctapipe.coordinates.CameraFrame

    Returns
    -------
    pandas.DataFrame

    """

    source_camera_pos = source.transform_to(camera_frame)
    data['source_x'] = source_camera_pos.x.value
    data['source_y'] = source_camera_pos.y.value

    return data


def get_sigma_time(
        data, theta2_on, theta2_off, theta2_cut=0.04, 
        norm_range=[0.2, 0.5], step_events=1):
    """
    Calculates time development of significance and background normalization

    Parameters
    ----------
    data: pandas.DataFrame with DL2/DL3 table
    theta2_on: pandas.DataFrame one column with theta^2 for the source position
    theta2_off: pandas.DataFrame one column with theta^2 for the OFF position
    theta2_cut: float
    norm_range: list
    step_events: int

    Returns
    -------
    list: sigma
    list: times
    list: obsid
    list: alphas

    """
    
    sigma = []
    alphas = []
    times = []
    obsid = []
    time = 0

    # number of off regions
    n_off = theta2_off.shape[1]

    mask_norm_on = (theta2_on > norm_range[0]) & (theta2_on < norm_range[1])
    mask_norm_off = (theta2_off > norm_range[0]) & (theta2_off < norm_range[1])
    mask_on = theta2_on <= theta2_cut
    mask_bkg = theta2_off <= theta2_cut

    for i in range(1, len(data['local_time']), step_events):

        mask_time = data['local_time'] <= data['local_time'][i]

        Non_norm = np.sum(mask_time & mask_norm_on)

        Noff_norm = 0
        for j in range(n_off):
            Noff_norm += np.sum(mask_norm_off[:, j] & mask_time)

        alpha = Non_norm / Noff_norm

        # this removes first event from each run, but I hope it is not a big deal
        if data['obs_id'][i]-data['obs_id'][i-1] == 0:
            N_ON = np.sum(mask_on & mask_time)

            N_OFF = 0
            for j in range(n_off):
                N_OFF += np.sum(mask_bkg[:, j] & mask_time)

            stat = WStatCountsStatistic(n_on=N_ON, n_off=N_OFF, alpha=alpha)
            sigma.append(stat.sqrt_ts)
            time += data['local_time'][i]-data['local_time'][i-1]
            times.append(time)
            obsid.append(data['obs_id'][i])
            alphas.append(alpha)

    # remove nans
    sigma = [0 if math.isnan(x) else x for x in sigma]
    alphas = [0 if math.isnan(x) else x for x in alphas]

    return sigma, times, obsid, alphas


def add_wobble_flag(data, horizon_frame=None, wobbles=[]):
    """
    Adds new column in the input table with wobble names based on the
    telescope pointing

    Parameters
    ----------
    data: astropy.table.Table

    horizon_frame: astropy.coordinates.AltAz
    wobbles: list of astropy.coordinates.SkyCoord

    Returns
    -------
    data: astropy.table.Table

    """

    pointing_coords = SkyCoord(alt = data['true_alt_tel'] * u.deg, 
                          az = data['true_az_tel'] * u.deg, 
                          frame=horizon_frame
                          )

    pointing_coords = pointing_coords.transform_to('icrs')
    data['wobble'] = 'W1'

    for wobble in wobbles:
        mask = pointing_coords.separation(wobble) < 0.1 * u.deg
        data['wobble'][mask] = wobble.info.name
    return data


def get_theta_off_stereo(data, n_off=1, on_region=None, wobbles=[]):
    """
    Calculates theta^2 for OFF regions regularly distributed around
    the pointing RA,DEC at the same offset as the ON region (for each 
    wobble). This is to be used for theta^2 in stereo, where the camera 
    frame is different for each telescope. The derivation of the rotation 
    angles is from:
    https://github.com/cta-observatory/magic-cta-pipe/blob/master/magicctapipe/utils/functions.py

    Parameters
    ----------
    data: astropy.table.Table
    n_off: int
    on_region: astropy.coordinates.SkyCoord
    wobbles: list of astropy.coordinates.SkyCoord

    Returns
    -------
    theta2_off: numpy.ndarray
        One column per OFF region
    off_radec: list of astropy.coordinates.SkyCoord
    
    """

    theta2_off = np.zeros([len(data), n_off])
    off_radec = []

    for wobble in wobbles:

        w_mask = data['wobble'] == wobble.info.name
        wobble_offset = on_region.separation(wobble)

        # Calculate the wobble rotation angle
        ra_diff = wobble.ra - on_region.ra

        numerator = np.sin(wobble.dec) * np.cos(on_region.dec)
        numerator -= np.cos(wobble.dec) * np.sin(on_region.dec) * np.cos(ra_diff)
        denominator = np.cos(wobble.dec) * np.sin(ra_diff)

        wobble_rotation = np.arctan2(numerator, denominator)
        wobble_rotation = Angle(wobble_rotation).wrap_at("360 deg")

        rotation_step = 360 / (n_off + 1)
        rotations_off = np.arange(0, 359, rotation_step) * u.deg

        rotations_off = rotations_off[rotations_off.to_value("deg") != 180]
        rotations_off += wobble_rotation

        off_coords = {}

        # Loop over every rotation angle
        for i_off, rotation in enumerate(rotations_off, start=1):
            skyoffset_frame = SkyOffsetFrame(origin=wobble, rotation=-rotation)

            # Calculate the OFF coordinate
            off_coord = SkyCoord(wobble_offset, "0 deg", frame=skyoffset_frame)
            off_coord = off_coord.transform_to("icrs")
            off_radec.append(off_coord)

            off_coords[i_off] = off_coord

        # Calculate the angular distance from the OFF coordinate
        for i_off, off_coord in off_coords.items():

            dl2_photon_list_off = data[w_mask].copy()
            event_coords = SkyCoord(data[w_mask]['reco_ra'] * u.deg, data[w_mask]['reco_dec'] * u.deg, frame='icrs')
            dl2_photon_list_off['theta'] = off_coord.separation(event_coords)

            dl2_photon_list_off['theta'] = dl2_photon_list_off['theta'].to(u.deg)

            theta2_off[w_mask, i_off-1] = np.array(dl2_photon_list_off['theta']**2)
        #print(theta2_off.shape)

    return theta2_off, off_radec


def get_theta_off(
        data, n_off=1, horizon_frame=None, 
        camera_frame=None, plot=False):
    """
    Calculates theta^2 for OFF regions regularly distributed around
    the center of the FoV at the same offset as the ON region. It
    also plots theta^2 for each off region.

    Parameters
    ----------
    data: astropy.table.Table
    n_off: int
    horizon_frame: astropy.coordinates.AltAz
    camera_frame: ctapipe.coordinates.CameraFrame
    plot: bool

    Returns
    -------
    theta2_off: numpy.ndarray
        One column per OFF region
    off_radec: list of astropy.coordinates.SkyCoord
    
    """

    true_source_position = [data['source_x'], data['source_y']]

    step = 2 * np.pi / (n_off + 1)
    angles = np.arange(0, 2 * np.pi, step)
    off_point = 0

    # ra dec of the off regions
    off_radec = []
    
    if plot:
        fig, ax = plt.subplots(1, len(angles)-1, figsize=(25,5))

    for i, angle in zip(range(len(angles)), angles):

        if angle > 0:
            x_off = np.cos(angle) * true_source_position[0] - np.sin(angle) * true_source_position[1]
            y_off = np.sin(angle) * true_source_position[0] + np.cos(angle) * true_source_position[1]

            off_source_position = [x_off, y_off]
            off_source_position[0].name = 'source_x'
            off_source_position[1].name = 'source_y'

            off_altaz = camera_to_altaz(
                np.array(off_source_position[0]) * u.m, 
                np.array(off_source_position[1]) * u.m, 
                horizon_frame=horizon_frame,
                camera_frame=camera_frame
            )
            
            wobbles = np.unique(np.array(data['wobble']))
            for w in wobbles:
                w_mask = data['wobble'] == w
                off_radec.append(off_altaz[w_mask][0].transform_to('icrs'))

            # theta for OFF coordinates, reflected background
            dl2_photon_list_off = get_theta(
                data.copy(), 
                zero_alt=off_altaz.alt.value, 
                zero_az=off_altaz.az.value
            )
            dl2_photon_list_off['theta'] = dl2_photon_list_off['theta'].to(u.deg)
            
            if plot:
                if len(angles) > 2:
                    axx = ax[i-1]
                else:
                    axx = ax
                h = axx.hist(data['theta']**2, bins=12, range=[0, 0.5], alpha=0.5, color='blue', label='ON')
                h2 = axx.hist(dl2_photon_list_off['theta']**2, bins=12, range=[0, 0.5], alpha=0.5, color='orange', label='OFF')
                axx.set_xlabel("$\\theta^{2} [deg^{2}]$")
                axx.legend()

            if off_point == 0:
                theta2_off = np.array(dl2_photon_list_off['theta']**2)
                theta2_off = np.expand_dims(theta2_off, axis=1)
            else:
                th2 = np.array(dl2_photon_list_off['theta']**2)
                theta2_off = np.concatenate((theta2_off, np.expand_dims(th2, axis=1)), axis=1)
            off_point += 1
    
    return theta2_off, off_radec


def camera_to_altaz(
        pos_x, pos_y, horizon_frame=None, 
        camera_frame=None):
    """
    Transforms coordinates in camera frame to horizon frame

    Parameters
    ----------
    pos_x: numpy.ndarray
    pos_y: numpy.ndarray
    horizon_frame: astropy.coordinates.AltAz
    camera_frame: ctapipe.coordinates.CameraFrame

    Returns
    -------
    horizon: astropy.coordinates.AltAz

    """

    camera_coord = SkyCoord(pos_x, pos_y, frame=camera_frame)
    horizon = camera_coord.transform_to(horizon_frame)

    return horizon


def get_theta2_from_dl3(dl3_path, good_obsids=None, target_coords=None, theta2_axis=None, n_off=5, norm_range=[0.5, 0.7]*u.deg, theta_cut=0.1*u.deg):

    data_store = DataStore.from_dir(dl3_path)
    theta2_off = np.zeros([len(theta2_axis.edges)-1, n_off])
    off_radec = []
    counts_all_on = []
    counts_all_all_off=[]

    event_counts = ThetaEventCounts()
    event_counts.n_off_regions = n_off

    sum_norm_on = 0
    sum_norm_off = 0
    N_on = 0
    N_off = 0
    t_elapsed = 0

    if good_obsids is not None:
        obsid_mask = [obsid in list(good_obsids.astype(int)) for obsid in data_store.obs_table["OBS_ID"]]
        good_obs_list = data_store.obs_table[obsid_mask]['OBS_ID']
        observations = data_store.get_observations(good_obs_list)
    else:
        observations = data_store.get_observations()

    for observation in observations:
        
        mask = data_store.obs_table['OBS_ID'] == observation.obs_id
        t_elapsed += data_store.obs_table[mask]['LIVETIME']

        # ON counts
        separation = target_coords.separation(observation.events.radec)
        
        N_on += sum(separation < theta_cut)
        
        counts_on, _ = np.histogram(separation ** 2, bins = theta2_axis.edges)
        counts_all_on.append(counts_on)
        
        norm_on = (separation > norm_range[0]) & (separation < norm_range[1])
        sum_norm_on += sum(norm_on)

        # OFF counts
        pos_angle = observation.pointing_radec.position_angle(target_coords)
        sep_angle = observation.pointing_radec.separation(target_coords)

        # Calculate the OFF counts from the wobble positions (OFF regions) provided
        rotation_step = 360 / (n_off + 1)
        rotations_off = np.arange(0, 359, rotation_step) * u.deg
        rotations_off = rotations_off[rotations_off.to_value("deg") != 0]
        rotations_off = pos_angle + rotations_off
        
        counts_all_off = []
        for i_off, rotation in enumerate(rotations_off, start=0):
            position_off = observation.pointing_radec.directional_offset_by(rotation, sep_angle)
            separation_off = position_off.separation(observation.events.radec)
            N_off += sum(separation_off < theta_cut)
            counts_off_wob, _ = np.histogram(separation_off ** 2, bins = theta2_axis.edges)
            norm_off = (separation_off > norm_range[0]) & (separation_off < norm_range[1])
            sum_norm_off += sum(norm_off)
            counts_all_off.append(counts_off_wob)

        counts_all_all_off.append(np.sum(np.array(counts_all_off), axis=0))
        
    alpha = sum_norm_on/sum_norm_off

    stat = WStatCountsStatistic(n_on=N_on, n_off=N_off, alpha=alpha)
    event_counts.significance_lima = stat.sqrt_ts
    event_counts.N_excess = N_on - alpha*N_off
    event_counts.t_elapsed = t_elapsed.to(u.h)[0]

    event_counts.N_on = N_on
    event_counts.N_off = N_off

    counts_on = np.sum(counts_all_on, axis=0)
    counts_off = np.sum(np.array(counts_all_all_off), axis=0)

    return counts_on, counts_off, alpha, event_counts


class ThetaEventCounts:

    def __init__(self):

        self.N_on = 0.
        self.N_off = 0.
        self.N_excess = 0.
        self.n_off_regions = 0.
        self.t_elapsed = 0.
        self.significance_lima = 0.