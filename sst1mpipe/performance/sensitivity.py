from ctapipe.io import read_table
from astropy.io.misc.hdf5 import write_table_hdf5
from astropy.table import Table, QTable
from astropy.coordinates import (
    angular_separation,
)
import numpy as np
import astropy.units as u
import tables
import ctaplot
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
import operator
from scipy import special
from scipy.interpolate import make_smoothing_spline

from pyirf.spectral import (
    calculate_event_weights,
    PowerLaw
)

from pyirf.cuts import (
    calculate_percentile_cut,
    evaluate_binned_cut
)

from gammapy.stats import WStatCountsStatistic

from sst1mpipe.io import (
    load_dl2_sst1m,
    check_outdir
)
from sst1mpipe.utils import (
    get_primary_type,
    correct_number_simulated_showers,
    mc_correct_shower_reuse,
    check_same_shower_fraction,
)

from .spectra import *


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

    particle_type = get_primary_type(input_file)
    simulated_event_info = QTable()

    for tel in ["tel_001", "stereo"]:
        try:
            params = read_table(input_file, "/dl2/event/telescope/parameters/" + tel)
        except Exception: pass
        try:
            params = read_table(input_file, "/dl1/event/telescope/parameters/" + tel)
        except Exception: pass
            

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
    simulated_event_info.viewcone = mc['max_viewcone_radius'][0] * mc['max_viewcone_radius'].unit
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

def relative_sensitivity(n_signal, n_background, alpha):
    """
    Calculates relative differential sensitivity in 
    units of Crab flux (usualy in a single energy bin).

    Parameters
    ----------
    n_signal: float
        Number of gammas in signal region
    n_background: float
        Number of gamma-like protons and 
        misreconstructed gammas in the background 
        regions.
    alpha: float
        Background normalization = 1/n_off_regions

    Returns
    ------- 
    n_excesses_5sigma: float
        Number of excess events to reach 5 sigma in
        given energy bin.
    sensitivity: float
        Percentage of Crab flux

    """

    stat = WStatCountsStatistic(
        n_on=n_signal + alpha * n_background,
        n_off=n_background,
        alpha=alpha
    )
    n_excesses_5sigma = stat.n_sig_matching_significance(5)
    sensitivity_5sigma_only = n_excesses_5sigma / (n_signal) * 100

    n_excesses_5sigma[n_excesses_5sigma < 10] = 10

    bkg_5percent = 0.05 * n_background * alpha
    if n_background > 0:
        n_excesses_5sigma[n_excesses_5sigma < bkg_5percent] = bkg_5percent[n_excesses_5sigma < bkg_5percent]

    sensitivity = n_excesses_5sigma / (n_signal) * 100  # percentage of Crab

    return n_excesses_5sigma, sensitivity_5sigma_only, sensitivity


def sensitivity_to_flux(
        sensitivity, energy, 
        target_spectrum=None):
    """
    Converts sensitivity in percent of give source 
    flux into differential sensitivity in flux units.

    Parameters
    ----------
    sensitivity: float
        Percentage of target spectrum flux
    energy: astropy.units.quantity.Quantity
        Center bin energy
    target_spectrum: pyirf.spectral.PowerLaw/pyirf.spectral.LogParabola/
    sst1mpipe.performance.spectra.PowerLawExpCutoff

    Returns
    ------- 
    astropy.units.quantity.Quantity

    """

    dFdE = target_spectrum(energy)
    sensitivity_flux = sensitivity / 100 * (dFdE * energy**2).to(u.TeV / (u.cm ** 2 * u.s))
    return sensitivity_flux


def get_time_to_detection(
        rate_signal=None, err_rate_signal=None, 
        rate_background=None, err_rate_background=None, 
        alpha=None, max_time=50 * u.hour):
    """
    Provides time needed to detection of given source 
    providing certain event rate on the background
    defined by event rate in the OFF regions.

    Parameters
    ----------
    rate_signal: astropy.units.quantity.Quantity
        Signal rate in the ON region (rate of gammas)
    err_rate_signal: astropy.units.quantity.Quantity
        Uncertainty on the signal rate
    rate_background: astropy.units.quantity.Quantity
        Background rate in the OFF regions
    err_rate_background: astropy.units.quantity.Quantity
        Uncertainty on the background rate
    alpha: float
        Background normalization = 1/n_off_regions
    max_time: astropy.units.quantity.Quantity
        Upper limit on time to detection above
        which the source is considered as 
        undetectable

    Returns
    ------- 
    astropy.units.quantity.Quantity
        Time to source detection
    astropy.units.quantity.Quantity
        Negative uncertainty on detection time
    astropy.units.quantity.Quantity
        Positive uncertainty on detection time

    """

    times = np.linspace(1, int(max_time.to_value(u.s)), int(max_time.to_value(u.s)/60.)) * u.s
    n_signal = rate_signal * times
    n_background = rate_background * times
    err_n_signal = err_rate_signal * times
    err_n_background = err_rate_background * times
    stat = WStatCountsStatistic(
        n_on=n_signal + alpha * n_background,
        n_off=n_background,
        alpha=alpha
    )
    stat_p = WStatCountsStatistic(
        n_on=n_signal + err_n_signal + alpha * (n_background - err_n_background),
        n_off=n_background - err_n_background,
        alpha=alpha
    )
    stat_m = WStatCountsStatistic(
        n_on=n_signal - err_n_signal + alpha * (n_background + err_n_background),
        n_off=n_background + err_n_background,
        alpha=alpha
    )
    mask = stat.sqrt_ts >= 5
    mask_p = stat_p.sqrt_ts >= 5
    mask_m = stat_m.sqrt_ts >= 5

    if sum(mask) == 0:
        logging.warning(f'5 sigma detection cannot be reached in {max_time.value} hours.')
        return max_time, 0 * u.hour, 0 * u.hour
    else:
        if sum(mask_p) == 0:
            err_minus = 0 * u.hour
        else:
            err_minus = min(times[mask]) - min(times[mask_p])
        if sum(mask_m) == 0:
            err_plus = 0 * u.hour
        else:
            err_plus = min(times[mask_m]) - min(times[mask])
        return min(times[mask]), err_minus, err_plus


def get_significance(n_signal=None, n_off=None, alpha=None):
    """
    Calculates statistical significance of the excess for 
    each pair of n_signal and n_off.

    Parameters
    ----------
    n_signal: numpy.ndarray
        Number of signal events per energy bin
    n_off: numpy.ndarray
        Number of background events per energy bin
    alpha: float
        Background normalization = 1/n_off_regions

    Returns
    ------- 
    numpy.ndarray
        Statistical significance of the excess
        per energy bin

    """

    stat = WStatCountsStatistic(
        n_on=n_signal + alpha * n_off,
        n_off=n_off,
        alpha=alpha
    )
    return stat.sqrt_ts


def check_spectrum(source):
    """
    Checks it the source spectrum is implemented 
    in sst1mpipe.performance.spectra

    Parameters
    ----------
    source: string
        E.g. \'CRAB_MAGIC_2020\'

    Returns
    ------- 

    """

    try:
        target_spectrum = globals()[source]
    except KeyError:
        logging.error('Desired spectrum is not implemented! \
        You can implement whatever you wish in \
        sst1mpipe.performance.spectra.')
        exit()


def get_gammaness_cuts(
        dl2_gamma, dl2_proton, config=None, 
        method='global', save_hdf=False, save_fig=False, 
        outdir=None, energy_bins=None, telescope=None, 
        gamma_off=False):
    """
    Provide gammaness cuts (masks) on input gamma and 
    proton tables. Different methods of cutting 
    may be selected.

    Parameters
    ----------
    dl2_gamma: astropy.table.Table
        DL2 MC point-like gamma table
    dl2_proton: astropy.table.Table
        DL2 MC diffuse proton table
    config: dict
    method: string
        Gammaness cut method to be applied \'global\',
        \'efficiency\', \'significance\'
    save_hdf: bool
        If True it stores energy dependent gammaness 
        cut table. It can be further used IRF production
        or in DL2->DL3 step in the data analysis
    save_fig: bool
        If True it stores a plot of energy dependent
        gammaness cuts
    outdir: string
    energy_bins: astropy.units.quantity.Quantity
    telescope: string
    gamma_off: bool
        If True it takes into account badly reconstructed 
        gammas which fall into the OFF region and 
        contribute to the background. It only works for 
        point-like gammas simulated with non-zero offset.
        So far it does not work for stereo sensitivity.

    Returns
    ------- 
    mask_gg:numpy.ndarray[bool]
        Gammaness mask for DL2 gamma table
    mask_gp: numpy.ndarray[bool]
        Gammaness mask for DL2 proton table
    """

    if method == 'global':

        gammaness_cut = config['analysis']['global_gammaness_cut']
        mask_gg = dl2_gamma['gammaness'] > gammaness_cut
        mask_gp = dl2_proton['gammaness'] > gammaness_cut
        logging.info('Global gammaness cut used: %f', gammaness_cut)

    elif method == 'efficiency':

        logging.info('Gammaness efficiency-based gammaness cut used.')
        min_events_bin = 100

        requested_gamma_efficiency = config['analysis']['gamma_efficiency']
        logging.info('Requested gamma efficiency %f', requested_gamma_efficiency)

        gammaness_cuts = calculate_percentile_cut(
                dl2_gamma["gammaness"],
                dl2_gamma["reco_energy"] * u.TeV,
                bins=energy_bins,
                min_value=0.1,
                max_value=0.95,
                fill_value=dl2_gamma["gammaness"].max(),
                percentile=100 * (1 - requested_gamma_efficiency),
                smoothing=None,
                min_events=min_events_bin,
            )

    elif method == 'significance':

        gammaness_cuts = calculate_gammaness_cuts_significance(dl2_gamma, dl2_proton, config=config, energy_bins=energy_bins, save_fig=save_fig, outdir=outdir, telescope=telescope, gamma_off=gamma_off)

    else:
        logging.error('Desired method of gammaness cut not implemented! Type \'sst1mpipe_mc_performance.py --help\' to see what is available.')
        exit()
       
    if save_hdf and (method != 'global'):
        logging.info(f'Differential sensitivity: saving energy dependent gammaness cut table for {telescope}..')

        outfile = outdir + '/gammaness_cuts_'+method+'_'+telescope+'.h5'
        write_table_hdf5(gammaness_cuts, outfile, path='gammaness_cuts', overwrite=True, append=True, serialize_meta=True)

    if save_fig and (method != 'global'):

        fig, ax = plt.subplots(1, 1, figsize=(8, 7))
        ax.errorbar(gammaness_cuts['center'], 
                    gammaness_cuts['cut'], 
                    xerr=(gammaness_cuts['center'] - gammaness_cuts['low'], gammaness_cuts['high'] - gammaness_cuts['center']),
                    fmt='o')

        energy_center = gammaness_cuts['center']
        ax.set_ylabel('gammaness cut')
        ax.set_xlabel(rf'$E_R$ [{energy_center.unit.to_string("latex")}]')
        ax.set_xscale('log')
        ax.grid(True, which='both')
        ax.set_xlim([min(gammaness_cuts['low'].value), max(gammaness_cuts['high'].value)])
        ax.set_ylim([0, 1.])
        fig.savefig(outdir + '/gammaness_cuts_'+method+'_'+telescope+'.png', dpi=200)

    if method != 'global':
        mask_gg = evaluate_binned_cut(
                dl2_gamma["gammaness"],
                dl2_gamma["reco_energy"] * u.TeV,
                gammaness_cuts,
                operator.ge,
            )

        mask_gp = evaluate_binned_cut(
                dl2_proton["gammaness"],
                dl2_proton["reco_energy"] * u.TeV,
                gammaness_cuts,
                operator.ge,
            )

    return mask_gg, mask_gp


def get_edep_theta_cuts(
        dl2_gamma, config=None, save_hdf=False, save_fig=False, 
        outdir=None, energy_bins=None, telescope=None, efficiency=0.68):
    """
    Calculates energy dependent theta2 cuts.

    Parameters
    ----------
    dl2_gamma: astropy.table.Table
        DL2 MC point-like gamma table
    config: dict
    save_hdf: bool
        If True it stores energy dependent gammaness 
        cut table. It can be further used IRF production
        or in DL2->DL3 step in the data analysis
    save_fig: bool
        If True it stores a plot of energy dependent
        gammaness cuts
    outdir: string
    energy_bins: astropy.units.quantity.Quantity
    telescope: string
    efficiency: float
       Requested fraction of events left after the cut in each energy bin. E.g. 0.68 gives angular resolution.

    For the moment works only in combination with global and efficiency gammaness cut

    Returns
    ------- 
    theta_cuts: astropy.table.Table

    """
    min_events_bin = 100
    requested_theta_efficiency = efficiency

    logging.info('Energy dependent theta2 cut used.')
    logging.info('Requested gamma efficiency: %f', efficiency)

    theta_cuts = calculate_percentile_cut(
            dl2_gamma["theta"].to(u.deg),
            dl2_gamma["reco_energy"] * u.TeV,
            bins=energy_bins,
            min_value=0.01 * u.deg,
            max_value=1.0 * u.deg,
            fill_value=0.5 * u.deg,
            percentile=100 * requested_theta_efficiency,
            smoothing=None,
            min_events=min_events_bin,
        )

    if save_hdf:

        outfile = outdir + '/theta_edep_cuts_'+telescope+'.h5'
        write_table_hdf5(theta_cuts, outfile, path='theta_cuts', overwrite=True, append=True, serialize_meta=True)

    if save_fig:

        fig, ax = plt.subplots(1, 1, figsize=(8, 7))
        ax.errorbar(theta_cuts['center'], 
                    theta_cuts['cut'], 
                    xerr=(theta_cuts['center'] - theta_cuts['low'], theta_cuts['high'] - theta_cuts['center']),
                    fmt='o')

        energy_center = theta_cuts['center']
        ax.set_ylabel('theta cut')
        ax.set_xlabel(rf'$E_R$ [{energy_center.unit.to_string("latex")}]')
        ax.set_xscale('log')
        ax.grid(True, which='both')
        ax.set_xlim([min(theta_cuts['low'].value), max(theta_cuts['high'].value)])
        ax.set_ylim([0, 1.])
        fig.savefig(outdir + '/theta_edep_cuts_'+telescope+'.png', dpi=200)

    theta2_cuts_bool = evaluate_binned_cut(
		dl2_gamma["theta"].to(u.deg),
		dl2_gamma["reco_energy"] * u.TeV,
		theta_cuts,
		operator.le,
		)

    return theta_cuts, theta_cuts_bool



def get_theta(dl2, zero_alt=None, zero_az=None):
    """
    Calculates theta^2 for ON region

    Parameters
    ----------
    dl2: astropy.table.Table
    zero_alt: astropy.table.column.Column
        Source altitude (column of floats in deg)
    zero_az: astropy.table.column.Column
        Source azimuth (column of floats in deg)

    Returns
    -------
    dl2: astropy.table.Table

    """

    # NOTE: Be carefull with the units! - must be revisited for proper unit conversion!
    if ~(zero_alt > np.pi).any():
        raise Exception("ERROR: true_alt, true_az are probably not in degrees! Check your DL2 file.")
    if ~(dl2['reco_alt'] > np.pi).any():
        raise Exception("ERROR: reco_alt, reco_az are probably not in degrees! Check your DL2 file.")

    dl2['theta'] = ctaplot.ana.angular_separation_altaz(
        dl2['reco_alt'] * np.pi/180. * u.rad, 
        dl2['reco_az'] * np.pi/180. * u.rad, 
        zero_alt * np.pi/180. * u.rad, 
        zero_az * np.pi/180. * u.rad
        )

    return dl2


def calculate_gammaness_cuts_significance(
        dl2_gamma, dl2_proton, config=None, 
        save_fig=False, outdir=None, energy_bins=None, 
        telescope=None, gamma_off=False):
    """
    Calculates energy dependent gammaness cut
    optimized to reach the best significance of 
    a source detection in each energy bin. This 
    cut estimator is biased towards detection of 
    a source with given spectrum and should be 
    used carefully. The DL2 input tables have to 
    already contain \'weight\' column to take 
    the source and background SED into account.

    Parameters
    ----------
    dl2_gamma: astropy.table.Table
        DL2 MC point-like gamma table
    dl2_proton: astropy.table.Table
        DL2 MC diffuse proton table
    config: dict
    save_fig: bool
        If True it stores a plot of energy dependent
        gammaness cuts
    outdir: string
    energy_bins: astropy.units.quantity.Quantity
    telescope: string
    gamma_off: bool
        If True it takes into account badly reconstructed 
        gammas which fall into the OFF region and 
        contribute to the background. It only works for 
        point-like gammas simulated with non-zero offset.
        So far it does not work for stereo sensitivity.

    Returns
    ------- 
    astropy.table.QTable

    """

    offset = angular_separation(dl2_gamma['true_az_tel'][0] * u.deg, dl2_gamma['true_alt_tel'][0] * u.deg, dl2_gamma['true_az'][0] * u.deg, dl2_gamma['true_alt'][0] * u.deg).to(u.deg)
    theta_cut = config['analysis']['global_theta_cut'] * u.deg
    theta_cut_p = 1.2 * u.deg
    area_ratio = theta_cut**2 / theta_cut_p**2
    N_off_regions = config['analysis']['off_regions']
    significance_fit = config['analysis']['significance_fit']

    gammaness_cuts = np.linspace(0.3, 0.95, 100)

    cut_table = QTable()
    cut_table["low"] = energy_bins[:-1]
    cut_table["high"] = energy_bins[1:]
    cut_table["n_gamma_events"] = 0.
    cut_table["n_proton_events"] = 0.
    cut_table["center"] = np.sqrt(energy_bins[:-1] * energy_bins[1:])
    cut_table["cut"] = dl2_gamma["gammaness"].max()

    for i in range(len(energy_bins)-1):
        mask_eg = (dl2_gamma['reco_energy'] * u.TeV > energy_bins[i]) & (dl2_gamma['reco_energy'] * u.TeV <= energy_bins[i+1])
        mask_ep = (dl2_proton['reco_energy'] * u.TeV > energy_bins[i]) & (dl2_proton['reco_energy'] * u.TeV <= energy_bins[i+1])

        n_signal_all = []
        n_off_all = []
        for gcut in gammaness_cuts:
            mask_gg = dl2_gamma['gammaness'] > gcut
            mask_gp = dl2_proton['gammaness'] > gcut

            mask = (dl2_gamma['theta'] < theta_cut) & mask_gg & mask_eg
            gammas_on = dl2_gamma[mask]
            n_signal_all.append(sum(gammas_on['weight']))

            mask = (dl2_proton['theta'] < theta_cut_p) & mask_gp & mask_ep
            protons_off = dl2_proton[mask]

            # get approximate number of gammas in off regions
            if offset > 0.01 * u.deg and gamma_off:
                if telescope == 'tel_001' or telescope == 'tel_002':
                    if sum(mask_eg & mask_gg) > 0:
                        n_gammas_off_observed = get_n_gammas_off(dl2_gamma[mask_eg & mask_gg], n_off=N_off_regions, theta_cut=theta_cut)
                    else:
                        n_gammas_off_observed = 0
            if gamma_off:
                n_off_all.append(sum(protons_off['weight']) * area_ratio.to_value(u.one) * N_off_regions + n_gammas_off_observed)
            else:
                n_off_all.append(sum(protons_off['weight']) * area_ratio.to_value(u.one) * N_off_regions)
        n_signal_all = np.array(n_signal_all)
        n_off_all = np.array(n_off_all)
    
        significance = get_significance(n_signal=n_signal_all, n_off=n_off_all, alpha=1/N_off_regions)
        if significance_fit:
            spl = make_smoothing_spline(gammaness_cuts, significance)
            y_plt = spl(gammaness_cuts)
            best_gcut = gammaness_cuts[np.argmax(y_plt)]
            logging.info('Getting energy dependent gammaness cut. Significance fit applied.')
        else:
            best_gcut = gammaness_cuts[np.argmax(significance)]
            logging.info('Getting energy dependent gammaness cut. Significance fit NOT applied')
            

        if save_fig and any(significance > 0):

            outdir_sens = outdir + '/g_cut_sens/'
            check_outdir(outdir_sens)
            emin = round(energy_bins[i].value, 2)
            emax = round(energy_bins[i+1].value, 2)
            fig, ax = plt.subplots(1, 1, figsize=(8, 7))
            ax.plot(gammaness_cuts, significance, color = 'b')
            if significance_fit:
                ax.plot(gammaness_cuts, y_plt, color = 'green')
            ax.axvline(x = best_gcut, color = 'r', label=f'best_gcut = {best_gcut}')
            ax.set_ylabel('significance')
            ax.set_xlabel('gammaness cut')
            ax.grid(True, which='both')
            ax.set_xlim([0, 1])
            ax.legend()
            fig.savefig(outdir_sens + '/significance_'+telescope+'_'+str(emin)+'_'+str(emax) + '.png', dpi=200)
            plt.close()
                
        cut_table["n_gamma_events"][i] = n_signal_all[np.argmax(significance)]
        cut_table["n_proton_events"][i] = n_off_all[np.argmax(significance)]
        if (n_signal_all[np.argmax(significance)] > 0) and (n_off_all[np.argmax(significance)] > 0):
            cut_table["cut"][i] = best_gcut

    return cut_table


def get_n_gammas_off(dl2_gamma, n_off=None, theta_cut=None):
    """
    Finds number of gammas which leak into the OFF 
    regions due to imperfect reconstruction

    Parameters
    ----------
    dl2_gamma: astropy.table.Table
        Point-gamma DL2 MC table
    n_off: int
        Number of OFF regions
    theta_cut: astropy.units.quantity.Quantity
        Cut on theta angle defining size of the 
        signal region and individual background
        regions.

    Returns
    -------
    int
        Number of gammas leaking in the OFF 
        regions

    """

    # this finds reflected regions wrt center of the FoV in camera coordinates
    # For stereo, it should be done based on reconstructed alt az or ra dec, and 
    # telescope pointing would be the average pointing of the two telescopes
    offset = angular_separation(
        dl2_gamma['true_az_tel'][0] * u.deg, 
        dl2_gamma['true_alt_tel'][0] * u.deg, 
        dl2_gamma['true_az'][0] * u.deg, 
        dl2_gamma['true_alt'][0] * u.deg
        ).to(u.deg)
    on_x = dl2_gamma['true_camera_x'][0]
    on_y = dl2_gamma['true_camera_y'][0]

    phi = 2*np.pi/(n_off+1)
    phi0 = np.arctan(on_y/on_x)

    focal = dl2_gamma['equivalent_focal_length'][0] * u.m
    theta_cut_m = np.tan(theta_cut) * focal

    N_gammas_off = 0
    for n in range(1, n_off+2):
        x = np.tan(np.cos(n*phi+phi0) * offset) * focal
        y = np.tan(np.sin(n*phi+phi0) * offset) * focal
        r = np.sqrt((x.to_value(u.m)-on_x)**2 + (y.to_value(u.m)-on_y)**2)
        # we need to skip the OFF region which overlaps with the ON region
        if r > 0.01:
            theta = np.sqrt((x.to_value(u.m)-dl2_gamma['reco_src_x'])**2 + (y.to_value(u.m)-dl2_gamma['reco_src_y'])**2)
            mask = (theta < theta_cut_m.to_value(u.m))
            gammas_off = dl2_gamma[mask]
            N_gammas_off += sum(gammas_off['weight'])
    return N_gammas_off


def sensitivity(
        input_file_gamma, input_file_proton, outdir=None, 
        config=None, telescope=None, save_fig=False, 
        save_hdf=False, gammaness_cuts=False, theta2_cuts=False, source_detection='', 
        energy_min=0.0, gamma_off=False):
    """
    Evaluates differential flux sensitivity on a point-like
    source. It also calculates optimal gammaness cuts
    and stores them in hdf table to be further used
    in IRF production and in DL2->DL3 on data. If spectrum
    of any point-like source is specified it provides 
    also estimated time to detection. Number of simulated
    events and re-used event fraction in each energy bin is
    listed to monitor sanity of pruduced results.

    Sensitivity in given energy bin is only computed if 
    N_simulated_gammas > 10 & N_simulated_protons > 10 
    in respective regions after the cuts. Otherwise the CTA 
    criteria on flux sensitivity in each energy bin are applied: 
        - significance >=5 sigma
        - excess > 5\% of background
        - N excess >= 10

    Parameters
    ----------
    input_file_gamma: string
        Path to DL2 testing MC point-like gammas
    input_file_proton: string
        Path to DL2 testing MC diffuse protons
    outdir: string
    config: dict
    telescope: string
    save_fig: bool
        If True some figures are stored
    save_hdf: bool
        If True hdf file with sensitivity
        table is stored
    gammaness_cuts: string
        Gammaness cut method to be applied \'global\',
        \'efficiency\', \'significance\'
    theta2_cuts: string
	Theta2 cut method to be applied \'global\', \'efficiency\'
    source_detection: string
        If not empty time to the source detection is 
        calculated. Spectrum must be specified in 
        performance.spectra
    energy_min: float
        Additional cut on minimum reconstructed energy
        in TeV. Only applied on estimated time to
        detect a source.
    gamma_off: bool
        If True it takes into account badly reconstructed 
        gammas which fall into the OFF region and 
        contribute to the background. It only works for 
        point-like gammas simulated with non-zero offset.
        So far it does not work for stereo sensitivity.

    Returns
    -------

    """

    logging.info('Estimating sensitivity for %s ...', telescope)

    obs_time = config['analysis']['observation_time_h'] * u.hour
    N_off_regions = config['analysis']['off_regions'] # assumed ratio between the size of the ON and OFF regions. 5 is fixed by CTA conventions for sensitivity calculation
    logging.info('N off regions %d', N_off_regions)
    logging.info('Assumed obstime %f h.', obs_time.value)

    energy_bins = np.logspace(
        config["analysis"]["log_energy_min_tev"], 
        config["analysis"]["log_energy_max_tev"],
        config["analysis"]["n_energy_bins"]
        ) * u.TeV

    dl2_gamma = load_dl2_sst1m(input_file_gamma, tel=telescope, config=config, table='astropy')
    dl2_proton = load_dl2_sst1m(input_file_proton, tel=telescope, config=config, table='astropy')

    mc_info_gamma = get_mc_info(input_file_gamma, config=config)
    mc_info_proton = get_mc_info(input_file_proton, config=config)

    # Traditionaly, the point gammas are weighted on Crab HEGRA Power-Law spectrum, not the MAGIC LogParabola,
    # but for the resulting flux sensitivity it doesn't matter.
    logging.info('CRAB_HEGRA spectrum is used for purposes of sensitivity estimation.')
    target_gamma_spectrum = CRAB_HEGRA # CRAB_HEGRA, CRAB_MAGIC_JHEAP2015
    #target_gamma_spectrum = CRAB_MAGIC_JHEAP2015
    dl2_gamma = get_weights(dl2_gamma, mc_info=mc_info_gamma, obs_time=obs_time, target_spectrum=target_gamma_spectrum)
    dl2_proton = get_weights(dl2_proton, mc_info=mc_info_proton, obs_time=obs_time, target_spectrum=DAMPE_P_He_SPECTRUM)

    # for point gammas, theta2 is calculated wrt the true simulated source position
    dl2_gamma = get_theta(dl2_gamma, zero_alt=dl2_gamma['true_alt'][0], zero_az=dl2_gamma['true_az'][0])
    # for diffuse protons theta2 is calculated wrt simulated point like gamma source position
    dl2_proton = get_theta(dl2_proton, zero_alt=dl2_gamma['true_alt'][0], zero_az=dl2_gamma['true_az'][0])

    # offset of point like gammas
    offset = angular_separation(dl2_gamma['true_az_tel'][0] * u.deg, dl2_gamma['true_alt_tel'][0] * u.deg, dl2_gamma['true_az'][0] * u.deg, dl2_gamma['true_alt'][0] * u.deg).to(u.deg)
    logging.info('Point gammas simulated with offset of: %f deg', offset.value)

    mask_gg, mask_gp = get_gammaness_cuts(dl2_gamma, dl2_proton, config=config, method=gammaness_cuts, save_hdf=save_hdf, save_fig=save_fig, outdir=outdir, energy_bins=energy_bins, telescope=telescope, gamma_off=gamma_off)

    if theta2_cuts == 'global':
	    theta_cut = config['analysis']['global_theta_cut'] * u.deg
	    logging.info('Global theta cut used: %f', theta_cut.value)
	    # Number of gammas in ON region after gammaness cut
	    mask = (dl2_gamma['theta'] < theta_cut) & mask_gg
    elif theta2_cuts == 'efficiency':
	    theta2_cuts_efficiency = config['analysis']['theta2_efficiency']
	    theta_cut, theta_cut_bool = get_edep_theta_cuts(dl2_gamma, config=config, save_hdf=save_hdf, save_fig=save_fig, outdir=outdir, energy_bins=energy_bins, telescope=telescope, efficiency=theta2_cuts_efficiency)
	    mask = theta_cut_bool & mask_gg
    else:
	    logging.error('Desired method of theta2 cut not implemented! Type \'sst1mpipe_mc_performance.py --help\' to see what is available')
	    exit()

    gammas_on = dl2_gamma[mask]
    logging.info('N of simulated gamma-like point gammas in ON region: %d', len(gammas_on))
    N_observed_g = sum(gammas_on['weight'])
    logging.info(f'Rate of observed gamma-like point gammas in ON region: {(N_observed_g/obs_time).to(1/u.s)}, Crab Spectrum')    

    # Number of protons in OFF region after gammaness cut
    # Here we take protons from a region centered at the same position as the ON region, 
    # but with larger radius to get better statistics of the background
    # Then to calculate the rate of remaining protons in ON region, we just use the ratio 
    # of areas of the two
    theta_cut_p = 1.2 * u.deg
    mask = (dl2_proton['theta'] < theta_cut_p) & mask_gp
    protons_off = dl2_proton[mask]
    logging.info(f'N of simulated gamma-like diffuse protons in OFF region (r={theta_cut_p}): {len(protons_off)}')

    # get approximate number of gammas in off regions
    if offset > 0.01 * u.deg and gamma_off:
        if telescope == 'tel_001' or telescope == 'tel_002':
            n_gammas_off_observed = get_n_gammas_off(dl2_gamma[mask_gg], n_off=N_off_regions, theta_cut=theta_cut)
            logging.info(f'Rate of observed gamma-like point gammas in all OFF regions: {(n_gammas_off_observed/obs_time).to(1/u.s)}, Crab Spectrum')   
        else:
            logging.info('Estimation of gammas in the OFF region is not yet supported.')

    if theta2_cuts == 'global':
	area_ratio = theta_cut**2 / theta_cut_p**2
	N_observed_p = sum(protons_off['weight']) * area_ratio
    elif theta2_cuts == 'efficiency':
	area_ratio = theta_cut['cut']**2 / theta_cut_p**2
	N_observed_p = 0
	for i in range(len(energy_bins)-1):
		maskp_theta = (protons_off['reco_energy'] * u.TeV > energy_bins[i]) & (protons_off['reco_energy'] * u.TeV <= energy_bins[i+1])
		N_observed_p += (sum(protons_off[maskp_theta]['weight']) * area_ratio[i])
    logging.info(f'Rate of observed gamma-like diffuse protons in ON region: {(N_observed_p/obs_time).to(1/u.s)}, p+He DAMPE')

    # Sanity check on the number of reused events in each energy bin and in ON/OFF regions after all cuts
    logging.info('GAMMA reused fractions:')
    fraction_g = check_same_shower_fraction(gammas_on.to_pandas(), energy_bins)
    logging.info('PROTON reused fractions:')
    fraction_p = check_same_shower_fraction(protons_off.to_pandas(), energy_bins)

    # Time to reach 5 sigma for all energies integrated
    if gamma_off:
        n_off = N_observed_p *  N_off_regions + n_gammas_off_observed
    else:
        n_off = N_observed_p *  N_off_regions

    time_5sig = get_time_to_detection(
        rate_signal=(N_observed_g/obs_time).to(1/u.s), 
        err_rate_signal=(np.sqrt(N_observed_g)/obs_time).to(1/u.s), 
        rate_background=(n_off/obs_time).to(1/u.s), 
        err_rate_background=(np.sqrt(n_off)/obs_time).to(1/u.s),
        alpha=1/N_off_regions
        )

    logging.info('Time to reach 5 sigma CRAB_HEGRA detection: %.2f +%.2f -%.2f hours', time_5sig[0].to_value(u.hour), time_5sig[2].to_value(u.hour), time_5sig[1].to_value(u.hour))

    # Get number of observed events in all energy bins
    sensitivity_all = np.zeros(len(energy_bins)-1)
    sensitivity_5sigma_only_all = np.zeros(len(energy_bins)-1)
    sensitivity_flux_all = np.zeros(len(energy_bins)-1)
    sensitivity_5sigma_only_flux_all = np.zeros(len(energy_bins)-1)
    sens_flux_err_minus = np.zeros(len(energy_bins)-1)
    sens_flux_err_plus = np.zeros(len(energy_bins)-1)
    energy = np.sqrt(energy_bins[:-1] * energy_bins[1:])

    logging.info("SENSITIVITY: [e min, e max] TeV, N sim gammas, N sim protons, N excess, N bkg, sensitivity [5sigma only] [\% Crab], sensitivity [\% Crab], sensitivity [TeV / (cm2 s)], sens err - [flux], sens err + [flux], fraction of G used N>1, fraction of P used N>1")
    for i in range(len(energy_bins)-1):
        maskg = (gammas_on['reco_energy'] * u.TeV > energy_bins[i]) & (gammas_on['reco_energy'] * u.TeV <= energy_bins[i+1])
        N_observed_g = sum(gammas_on[maskg]['weight'])
        N_simulated_g = len(gammas_on[maskg])

        maskp = (protons_off['reco_energy'] * u.TeV > energy_bins[i]) & (protons_off['reco_energy'] * u.TeV <= energy_bins[i+1])
	
	if theta_cuts == 'global':
        	N_observed_p = sum(protons_off[maskp]['weight']) * area_ratio
        elif theta2_cuts == 'efficiency':
		N_observed_p = sum(protons_off[maskp]['weight']) * area_ratio[i]
	
	N_simulated_p = len(protons_off[maskp])

        # get approximate number of gammas in off regions
        if offset > 0.01 * u.deg and gamma_off:
            if telescope == 'tel_001' or telescope == 'tel_002':
                mask_e = (dl2_gamma['reco_energy'] * u.TeV > energy_bins[i]) & (dl2_gamma['reco_energy'] * u.TeV <= energy_bins[i+1])
                if sum(mask_gg & mask_e) > 0:
                    n_gammas_off_observed = get_n_gammas_off(dl2_gamma[mask_gg & mask_e], n_off=N_off_regions, theta_cut=theta_cut)
                else:
                    n_gammas_off_observed = 0

        n_signal = N_observed_g
        if gamma_off:
            n_off = n_gammas_off_observed + N_observed_p * N_off_regions
        else:
            n_off = N_observed_p * N_off_regions
        alpha = 1/N_off_regions

        if (N_simulated_g > 10) & (N_simulated_p > 10) & (n_signal >= 0) & (n_off >= 0):
            # relative sensitivity in % or Crab (5sigma, excess > 5\% of background)
            _, sensitivity_5sigma_only, sensitivity = relative_sensitivity(n_signal, n_off, alpha)

            # Convert sensitivity to flux
            sensitivity_flux = sensitivity_to_flux(sensitivity, energy[i], target_spectrum=target_gamma_spectrum)
            sensitivity_5sigma_only_flux = sensitivity_to_flux(sensitivity_5sigma_only, energy[i], target_spectrum=target_gamma_spectrum)

            # Uncertainty on sensitivity taking into account uncertainty in signal and backgound
            # We cannot calculate Poisson uncertainty from re-weighted events, because then it would be source dependent.
            # According to Vladimir, we should calculate it from simulated events
            err_n_signal = np.sqrt(N_simulated_g) * n_signal / N_simulated_g
            err_n_off = np.sqrt(N_simulated_p) * n_off / N_simulated_p

            if n_off < err_n_off:
                _, _, sens_better = relative_sensitivity(n_signal + err_n_signal, 0, alpha)
            else:
                _, _, sens_better = relative_sensitivity(n_signal + err_n_signal, n_off - err_n_off, alpha)
            _, _, sens_worse = relative_sensitivity(n_signal - err_n_signal, n_off + err_n_off, alpha)
            sens_flux_better = sensitivity_to_flux(sens_better, energy[i], target_spectrum=target_gamma_spectrum)
            sens_flux_worse = sensitivity_to_flux(sens_worse, energy[i], target_spectrum=target_gamma_spectrum)

        else:
            sensitivity = np.nan
            sensitivity_5sigma_only = np.nan
            sensitivity_flux = np.nan * u.TeV / (u.cm ** 2 * u.s)
            sensitivity_5sigma_only_flux = np.nan * u.TeV / (u.cm ** 2 * u.s)
            sens_flux_better = np.nan * u.TeV / (u.cm ** 2 * u.s)
            sens_flux_worse = np.nan * u.TeV / (u.cm ** 2 * u.s)

        sensitivity_all[i] = sensitivity
        sensitivity_5sigma_only_all[i] = sensitivity_5sigma_only
        sensitivity_flux_all[i] = sensitivity_flux.value
        sensitivity_5sigma_only_flux_all[i] = sensitivity_5sigma_only_flux.value
        sens_flux_err_minus[i] = sensitivity_flux.value-sens_flux_better.value
        sens_flux_err_plus[i] = sens_flux_worse.value-sensitivity_flux.value
        if sens_flux_err_minus[i] < 0:
            sens_flux_err_minus[i] = 0.0
            logging.warning(f'[{energy_bins[i]:.2f}, {energy_bins[i+1]:.2f}], flux_sensitivity_err_minus cannot be estimated!')
        if sens_flux_err_plus[i] < 0:
            sens_flux_err_plus[i] = 10*sensitivity_flux.value
            logging.warning(f'[{energy_bins[i]:.2f}, {energy_bins[i+1]:.2f}], flux_sensitivity_err_plus cannot be estimated!')
        # We print also simulated (not weighted) number of protons and gammas after all cuts to see how good is the statistics 
        # on which the sensitivity is based
        logging.info(f'[{energy_bins[i].to_value(u.TeV):6.1f}, {energy_bins[i+1].to_value(u.TeV):6.1f}], {N_simulated_g:4d}, {N_simulated_p:4d}, {n_signal:5.1f}, {n_off:5.1f}, {sensitivity_5sigma_only:5.1f}, {sensitivity:5.1f}, {sensitivity_flux.value:3.1E}, {sens_flux_err_minus[i]:3.1E}, {sens_flux_err_plus[i]:3.1E}, {fraction_g[i, 2]:3.1E}, {fraction_p[i, 2]:3.1E}')
    
    # NOTE: This is uggly. We get rid of the units in order to merge everything in an array, and later adding the units after conversion
    # in the Astropy Table. 
    data = np.array(np.column_stack((energy.value, energy_bins[:-1].value, energy_bins[1:].value, sensitivity_all, sensitivity_flux_all, sens_flux_err_minus, sens_flux_err_plus)))

    names = ['energy', 'energy_min', 'energy_max', 'relative_sensitivity', 'flux_sensitivity', 'flux_sensitivity_err_minus', 'flux_sensitivity_err_plus']
    sens = Table(data=data, names=names)
    sens['energy'].unit = energy.unit
    sens['energy_min'].unit = energy_bins.unit
    sens['energy_max'].unit = energy_bins.unit
    sens['flux_sensitivity'].unit = sensitivity_flux.unit

    if save_fig:

        logging.info(f'Differential sensitivity: saving images for {telescope}..')
        fig = plt.figure(figsize=(6, 5))
        
        plt.ylabel(rf'Flux sensitivity [{sensitivity_flux.unit.to_string("latex")}]')
        plt.xlabel(rf'$E_R$ [{energy.unit.to_string("latex")}]')
        plt.xscale('log')
        plt.yscale('log')

        plt.errorbar(sens['energy'], sens['flux_sensitivity'],
                    xerr=(sens['energy'] - sens['energy_min'],sens['energy_max'] - sens['energy']),
                    yerr=(sens['flux_sensitivity_err_minus'],sens['flux_sensitivity_err_plus']),
                    fmt='o',
                    label='Diff. sensitivity of ' + telescope
                    )
        plt.errorbar(energy, sensitivity_5sigma_only_flux_all,
                    xerr=(energy - energy_bins[:-1], energy_bins[1:] - energy),
                    fmt='.',
                    alpha=0.5,
                    label='5sigma condition only'
                    )
        # Crab
        energy_smooth = np.logspace(-1, 3, 200) * u.TeV
        crab_flux = CRAB_MAGIC_JHEAP2015(energy_smooth)
        plt.plot(energy_smooth, (crab_flux * energy_smooth**2).to(u.TeV / (u.cm ** 2 * u.s)), color='grey', label="Crab (JHEAP2015)")

        plt.xlim([10**-1, 5*10**2])
        plt.ylim([10**-12, 10**-10])
        plt.grid(True, which='both')
        plt.legend()

        fig.tight_layout()
        fig.savefig(outdir + '/flux_sensitivity_'+telescope+'.png', dpi=200)

    if save_hdf:
        logging.info(f'Differential sensitivity: saving tables for {telescope}..')

        outfile = outdir + '/flux_sensitivity_'+telescope+'.h5'
        write_table_hdf5(sens, outfile, path='sensitivity', overwrite=True, append=True, serialize_meta=True)

    # estimation of time needed for detection of a given source
    if len(source_detection) > 0:
        
        dl2_gamma.remove_column('weight')
        dl2_proton.remove_column('weight')
        check_spectrum(source_detection)
        source_time_to_detection(
            dl2_gamma, 
            dl2_proton, 
            mc_info_gamma=mc_info_gamma, 
            mc_info_proton=mc_info_proton, 
            config=config, 
            telescope=telescope, 
            gammaness_cuts=gammaness_cuts, 
            theta2_cuts=theta2_cuts,
            source=source_detection, 
            energy_min=energy_min)


def source_time_to_detection(
        dl2_gamma, dl2_proton, mc_info_gamma=None, 
        mc_info_proton=None, config=None, telescope=None, 
        gammaness_cuts=None, theta2_cuts=None,
	source=None, energy_min=0.0):
    """
    Provides time needed to detection of given source above
    certain energy.

    Here we basicaly recalculate sensitivity from DL2 MC tables.
    In the future, this should be calculated with the use of 
    previously calculated IRFs.

    Parameters
    ----------
    dl2_gamma: astropy.table.Table
        Point-gamma DL2 MC table
    dl2_proton: astropy.table.Table
        Diffuse proton DL2 MC table
    mc_info_gamma: astropy.table.QTable
        Extracted info about MC production
    mc_info_proton: astropy.table.QTable
        Extracted info about MC production
    config: dict
    telescope: string
    gammaness_cuts: string
        Gammaness cut method to be applied \'global\',
        \'efficiency\', \'significance\'
    theta2_cuts: string
	Theta2 cut method to be applied \'global\', \'efficiency\'
    source: string
        Source spectrum defined in 
        sst1mpipe.performance.spectra
    energy_min: float
        Additional cut on minimum reconstructed energy
        in TeV. Only applied on estimated time to
        detect a source.

    Returns
    -------

    """

    logging.info('Estimating time to detection of %s for %s ...', source, telescope)

    obs_time = config['analysis']['observation_time_h'] * u.hour
    N_off_regions = config['analysis']['off_regions'] # assumed ratio between the size of the ON and OFF regions. 5 is fixed by CTA conventions for sensitivity calculation
    logging.info('N off regions %d', N_off_regions)

    energy_bins = np.logspace(
        config["analysis"]["log_energy_min_tev"], 
        config["analysis"]["log_energy_max_tev"],
        config["analysis"]["n_energy_bins"]
        ) * u.TeV

    dl2_gamma = dl2_gamma[dl2_gamma['reco_energy'] >= energy_min]
    dl2_proton = dl2_proton[dl2_proton['reco_energy'] >= energy_min]
    logging.info('Cut on minumum reco energy > %f TeV applied.', energy_min)

    # here we magicaly convert the string into the name of a function
    # https://bobbyhadz.com/blog/python-call-function-by-string-name
    target_spectrum = globals()[source]
    dl2_gamma = get_weights(dl2_gamma, mc_info=mc_info_gamma, obs_time=obs_time, target_spectrum=target_spectrum)
    dl2_proton = get_weights(dl2_proton, mc_info=mc_info_proton, obs_time=obs_time, target_spectrum=DAMPE_P_He_SPECTRUM)

    # for point gammas, theta2 is calculated wrt the true simulated source position
    dl2_gamma = get_theta(dl2_gamma, zero_alt=dl2_gamma['true_alt'][0], zero_az=dl2_gamma['true_az'][0])
    # for diffuse protons theta2 is calculated wrt simulated point like gamma source position
    dl2_proton = get_theta(dl2_proton, zero_alt=dl2_gamma['true_alt'][0], zero_az=dl2_gamma['true_az'][0])

    mask_gg, mask_gp = get_gammaness_cuts(dl2_gamma, dl2_proton, config=config, method=gammaness_cuts, energy_bins=energy_bins, telescope=telescope)

    if theta2_cuts == 'global':
        theta_cut = config['analysis']['global_theta_cut'] * u.deg
        logging.info('Global theta cut used: %f', theta_cut.value)
        # Number of gammas in ON region after gammaness cut
        mask = (dl2_gamma['theta'] < theta_cut) & mask_gg
    elif theta2_cuts == 'efficiency':
	    theta2_cuts_efficiency = config['analysis']['theta2_efficiency']
	    theta_cut, theta_cut_bool = get_edep_theta_cuts(dl2_gamma, config=config, energy_bins=energy_bins, telescope=telescope, efficiency=theta2_cuts_efficiency)
	    mask = theta_cut_bool & mask_gg
    else:
	    logging.error('Desired method of theta2 cut not implemented! Type \'sst1mpipe_mc_performance.py --help\' to see what is available.')
	    exit()


    gammas_on = dl2_gamma[mask]
    logging.info('N of simulated gamma-like point gammas in ON region: %d', len(gammas_on))
    N_observed_g = sum(gammas_on['weight'])
    logging.info(f'Rate of observed gamma-like point gammas in ON region: {(N_observed_g/obs_time).to(1/u.s)}, {source} Spectrum')    

    # Number of protons in OFF region after gammaness cut
    # Here we take protons from a region centered at the same position as the ON region, 
    # but with larger radius to get better statistics of the background
    # Then to calculate the rate of remaining protons in ON region, we just use the ratio 
    # of areas of the two
    theta_cut_p = 1.2 * u.deg
    mask = (dl2_proton['theta'] < theta_cut_p) & mask_gp
    protons_off = dl2_proton[mask]
    logging.info(f'N of simulated gamma-like diffuse protons in OFF region (r={theta_cut_p}): {len(protons_off)}')


    if theta2_cuts == 'global':
	    area_ratio = theta_cut**2 / theta_cut_p**2
	    N_observed_p = sum(protons_off['weight']) * area_ratio
    elif theta2_cuts == 'efficiency':
	    area_ratio = theta_cut['cut']**2 / theta_cut_p**2
	    N_observed_p = 0
	for i in range(len(energy_bins)-1):
		maskp_theta = (protons_off['reco_energy'] * u.TeV > energy_bins[i]) & (protons_off['reco_energy'] * u.TeV <= energy_bins[i+1])
		N_observed_p += (sum(protons_off[maskp_theta]['weight']) * area_ratio[i])

    logging.info(f'Rate of observed gamma-like diffuse protons in ON region: {(N_observed_p/obs_time).to(1/u.s)}, p+He DAMPE')

    # Time to reach 5 sigma for all energies integrated
    time_5sig = get_time_to_detection(
        rate_signal=(N_observed_g/obs_time).to(1/u.s), 
        err_rate_signal=(np.sqrt(N_observed_g)/obs_time).to(1/u.s), 
        rate_background=(N_observed_p/obs_time).to(1/u.s) *  N_off_regions, 
        err_rate_background=(np.sqrt(N_observed_p * N_off_regions)/obs_time).to(1/u.s),
        alpha=1/N_off_regions,
        max_time=500 * u.hour
        )
    logging.info('Time to reach 5 sigma %s detection: %.2f +%.2f -%.2f hours', source, time_5sig[0].to_value(u.hour), time_5sig[2].to_value(u.hour), time_5sig[1].to_value(u.hour))
