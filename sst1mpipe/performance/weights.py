"""
Weights of the simulated events (MC information of the files, spectral reweighting)
"""
import logging

import astropy.units as u
from astropy.table import QTable
from ctapipe.io import read_table
from pyirf.spectral import PowerLaw, calculate_event_weights

from sst1mpipe.utils import correct_number_simulated_showers, mc_correct_shower_reuse

__all__ = ["get_mc_info", "get_weights"]


def get_mc_info(input_file, config=None):
    """
    Extracts MC simulation setup (energy range, N thrown
    showers, etc.) from the input DL1/DL2 MC files and
    applies some non trivial corrections which are necessary
    due to a small bug in CORSIKA or when DL1 resulting from
    several different MC productions are combined.

    Parameters
    ----------
    input_file: string
        Path to DL1/DL2 MC file
    config: dict

    Returns
    -------
    astropy.table.QTable

    """

    logging.info('Getting MC info for %s', input_file)

    mc = read_table(input_file, "/configuration/simulation/run")
    histograms = read_table(input_file, "/simulation/service/shower_distribution")

    # Interestingly, if we use CSCAT>20 the number stored in mc['shower_reuse'] is still 20, even though the
    # total number of shower stored in the simtel file is correct (Vladimir checked), and also the total
    # number of thrown events stored in histograms seems to be fine (Jakub checked). Here we use the total
    # number of events stored in histograms to correct for shower_reuse. NOTE: should be implemented
    # at the r1->dl1 level in the future!
    mc = mc_correct_shower_reuse(mc, histograms)

    simulated_event_info = QTable()

    for tel in ["tel_001", "stereo"]:
        try:
            params = read_table(input_file, "/dl2/event/telescope/parameters/" + tel)
        except Exception:
            pass
        try:
            params = read_table(input_file, "/dl1/event/telescope/parameters/" + tel)
        except Exception:
            pass


    if 'min_true_energy_cut' in params.keys():

        min_energy_cfg = params['min_true_energy_cut'][0]

        if (min_energy_cfg > mc['energy_range_min']).any():
            logging.warning('There was a cut on energy applied on this DL2 file, which doesn''t correspond with \'energy_range_min\' in MC DL2 tab).\nTHE TOTAL NUMBER OF SIMULATED EVENTS WILL BE RECALCULATED FOR SPECTRAL REWEIGHTING.')
            simulated_event_info.energy_min = min_energy_cfg * u.TeV

            # Correcting the number of totaly simulated showers after the cuts on true energy in DL1 and DL2
            simulated_event_info = correct_number_simulated_showers(simulated_event_info, mc_table=mc, histogram_table=histograms, e_min_cut=min_energy_cfg)
        else:
            simulated_event_info.energy_min = mc['energy_range_min'][0] * mc['energy_range_min'].unit
            simulated_event_info.n_showers = (mc['n_showers'] * mc['shower_reuse']).sum()
    else:
        simulated_event_info.energy_min = mc['energy_range_min'][0] * mc['energy_range_min'].unit
        simulated_event_info.n_showers = (mc['n_showers'] * mc['shower_reuse']).sum()

    simulated_event_info.energy_max = mc['energy_range_max'][0] * mc['energy_range_max'].unit
    simulated_event_info.spectral_index = mc['spectral_index'][0]
    simulated_event_info.viewcone_min = mc['min_viewcone_radius'][0] * mc['min_viewcone_radius'].unit
    simulated_event_info.viewcone_max = mc['max_viewcone_radius'][0] * mc['max_viewcone_radius'].unit
    simulated_event_info.min_impact = mc['min_scatter_range'][0] * mc['min_scatter_range'].unit
    simulated_event_info.max_impact = mc['max_scatter_range'][0] * mc['max_scatter_range'].unit

    return simulated_event_info


def get_weights(
        data, mc_info=None, obs_time=None,
        target_spectrum=None):
    """
    Calculate event weights to MC event distribution for
    any source. Weights are stored as extra column in the
    data table.

    Parameters
    ----------
    data: astropy.table.Table
        DL1/DL2 MC table
    mc_info: astropy.table.QTable
        Table with information about simulated distribution
        of MC events
    obs_time: astropy.units.quantity.Quantity
        Target observing time
    target_spectrum: pyirf.spectral.PowerLaw/pyirf.spectral.LogParabola/
    sst1mpipe.performance.spectra.PowerLawExpCutoff

    Returns
    -------
    astropy.table.Table

    """

    simulated_spectrum = PowerLaw.from_simulation(mc_info, obs_time)

    if (data["true_energy"] * u.TeV).unit is not u.TeV:
        raise Exception('ERROR: DL2 table has apparently stored units now! Remove unit addition for true_energy in get_weights()')

    data["weight"] = calculate_event_weights(
            data["true_energy"] * u.TeV,
            target_spectrum,
            simulated_spectrum
            )

    return data
