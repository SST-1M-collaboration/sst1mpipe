from ctapipe.io import read_table
from astropy.io.misc.hdf5 import write_table_hdf5
from astropy.table import Table, QTable
import numpy as np
import pandas as pd
import astropy.units as u
import tables
import ctaplot
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn import metrics
import logging

import shutil

from pathlib import Path

from astropy.io import fits
from astropy.coordinates import AltAz

import sst1mpipe

from sst1mpipe.io import(
    load_dl2_sst1m,
    check_outdir,
    load_config
)
from sst1mpipe.utils import (
    check_same_shower_fraction,
    get_avg_pointing,
)

from .spectra import DAMPE_P_He_SPECTRUM, CRAB_HEGRA
from .sensitivity import (
    get_mc_info,
    get_weights,
    get_gammaness_cuts,
    get_theta
)

from pyirf.simulations import SimulatedEventsInfo

from pyirf.irf import (effective_area, 
                       effective_area_per_energy,
                       effective_area_per_energy_and_fov, 
                       psf_table,
                       background_2d,
                       energy_dispersion
                      )
from pyirf.io import (
    create_aeff2d_hdu,
    create_psf_table_hdu,
    create_energy_dispersion_hdu,
    create_rad_max_hdu,
    create_background_2d_hdu,
)
from astropy.io.misc.hdf5 import read_table_hdf5
from pyirf.cuts import (
    evaluate_binned_cut
)
import operator


def evaluate_performance(
        gamma_file=None, proton_file=None, 
        outdir=None, config=None, telescope=None, 
        save_fig=False, save_hdf=False, 
        gammaness_cuts=None):
    """
    Random Forrest performance evaluation on DL2 MC testing
    point-like gammas and diffuse protons. It provides 
    energy and angular resolution and if proton file is 
    provided (optional) ROC curve is also evaluated.

    Except for global gammaness cut option, resolutions
    are calculated for events after optimal energy dependend
    cut (the same that is used for sensitivity curve),
    which corresponds to the real observation of 
    a point-like source.

    Parameters
    ----------
    gamma_file: string
        Path to DL2 MC point-like gammas
    proton_file: string
        Path to DL2 MC diffuse protons
    outdir: string
    config: dict
    telescope: string
    save_fig: bool
        If True some figures are stored
    save_hdf: bool
        If True hdf files with performance
        tables are stored
    gammaness_cuts: string
        Gammaness cut method to be applied \'global\',
        \'efficiency\', \'significance\'

    Returns
    -------

    """

    logging.info('Evaluating performance for %s', telescope)

    energy_bins = np.logspace(
    config["analysis"]["log_energy_min_tev"], 
    config["analysis"]["log_energy_max_tev"],
    config["analysis"]["n_energy_bins"]
    ) * u.TeV

    dl2_gamma = load_dl2_sst1m(gamma_file, tel=telescope, config=config, table='pandas')

    if proton_file is not None:
        dl2_proton = load_dl2_sst1m(proton_file, tel=telescope, config=config, table='pandas')

        # Mixing gammas and protons
        gh_testing_dataset = pd.concat([dl2_gamma, dl2_proton], ignore_index=True)
        gh_testing_dataset = gh_testing_dataset.sample(frac=1).reset_index(drop=True)

    # select only events that pass a global gammaness cut
    if (dl2_gamma['gammaness'] == 0).all():
        logging.info('Gammaness of all events in the input file is equal to 0. Gamma-hadron classifier probably was not applied in DL1->DL2 step and thus no gammaness cut is applied here!')
    else:
        if gammaness_cuts != 'global':
            logging.info('Energy and angular resolution are calculated only for gammas which passed ENERGY DEPENDENT GAMMANESS CUT. Optimized on %s', gammaness_cuts)

            dl2_gamma_table = Table.from_pandas(dl2_gamma)
            dl2_proton_table = Table.from_pandas(dl2_proton) 

            if gammaness_cuts == 'significance':

                # for point gammas, theta2 is calculated wrt the true simulated source position
                dl2_gamma_table = get_theta(dl2_gamma_table, zero_alt=dl2_gamma_table['true_alt'][0], zero_az=dl2_gamma_table['true_az'][0])
                # for diffuse protons theta2 is calculated wrt simulated point like gamma source position
                dl2_proton_table = get_theta(dl2_proton_table, zero_alt=dl2_proton_table['true_alt'][0], zero_az=dl2_proton_table['true_az'][0])
                
                mc_info_gamma = get_mc_info(gamma_file, config=config)
                mc_info_proton = get_mc_info(proton_file, config=config)

                obs_time = 50 * u.hour
                dl2_gamma_table = get_weights(dl2_gamma_table, mc_info=mc_info_gamma, obs_time=obs_time, target_spectrum=CRAB_HEGRA)
                dl2_proton_table = get_weights(dl2_proton_table, mc_info=mc_info_proton, obs_time=obs_time, target_spectrum=DAMPE_P_He_SPECTRUM)

            mask, _ = get_gammaness_cuts(dl2_gamma_table, dl2_proton_table, config=config, method=gammaness_cuts, telescope=telescope)
            dl2_gamma = dl2_gamma[mask].copy()
        else:
            logging.info('Energy and angular resolution are calculated only for gammas which passed GLOBAL GAMMANESS > %f cut.', config["analysis"]["global_gammaness_cut"])
            mask = dl2_gamma["gammaness"] > config["analysis"]["global_gammaness_cut"]
            dl2_gamma = dl2_gamma[mask].copy()

    # NOTE: I am not sure if these should be binned in true (as we do in the lst performance paper) 
    # of reconstructed energy (as ctaplot seems to be suggesting in some functions). For now there is 
    # a switch. Default is true energy binning.
    energy_resolution(dl2_gamma, e_bins=energy_bins, outdir=outdir, telescope=telescope, save_fig=save_fig, save_hdf=save_hdf, x_axis_true_energy=True)

    # In mono, here we apply cut on wrongly reconstructed disp sign (following the LST performance paper)
    # Events with wrongly reco sign are mostly on low energies and they create
    # a donut shape distorsion of the PSF. The radius of that donut is about 1.5 deg
    # so it does not affect our ability to distinguish two point-like sources.
    if telescope == 'stereo':
        logging.info('Telescope ' + telescope + '. Cut on disp sign is NOT APPLIED for evaluation of angular resolution.')
        angular_resolution(
            dl2_gamma, e_bins=energy_bins, outdir=outdir, 
            telescope=telescope, save_fig=save_fig, save_hdf=save_hdf, 
            x_axis_true_energy=True, axes_sky=True
            )
    else:
        logging.info('Telescope ' + telescope + '. Cut on disp sign IS APPLIED for evaluation of angular resolution.')
        mask_disp_sign = dl2_gamma['reco_disp_sign'] == dl2_gamma['disp_sign']
        angular_resolution(
            dl2_gamma[mask_disp_sign], e_bins=energy_bins, outdir=outdir, 
            telescope=telescope, save_fig=save_fig, save_hdf=save_hdf, 
            x_axis_true_energy=True, axes_sky=True
            )
    
    if proton_file is not None:
        roc_curve(gh_testing_dataset, e_bins=energy_bins, outdir=outdir, telescope=telescope, save_fig=save_fig, save_hdf=save_hdf)
    else:
        logging.warning('DL2 proton file not specified. Performance of RF classifier (ROC curve) cannot be estimated.')


def energy_resolution_per_energy(
        true_energy, reco_energy, percentile=68.27, 
        confidence_level=0.95, bias_correction=False, 
        bins=None, x_axis_true_energy=True):
    """
    Modified functions from ctaplot, which doesnt 
    seem to work properly. Evaluates energy resolution 
    per bin in energy

    Parameters
    ----------
    true_energy: astropy.units.quantity.Quantity
    reco_energy: astropy.units.quantity.Quantity
    percentile: float
    confidence_level: float
    bias_correction: bool
    bins: astropy.units.quantity.Quantity
        Energy bin edges
    x_axis_true_energy: bool
        If True energy axis is binned in true
        energy. Reconstructed energy is used otherwise

    Returns
    -------
    bins: astropy.units.quantity.Quantity
    res_e: numpy.ndarray

    """

    res_e = []

    for i in range(len(bins)-1):
        if x_axis_true_energy:
            mask = (true_energy > bins[i]) & (true_energy < bins[i + 1])
        else:
            mask = (reco_energy > bins[i]) & (reco_energy < bins[i + 1])

        res_e.append(ctaplot.ana.energy_resolution(true_energy[mask], reco_energy[mask],
                                       percentile=percentile,
                                       confidence_level=confidence_level,
                                       bias_correction=bias_correction))

    res_e = np.array(res_e, dtype=float)
    # replace zeros with Nans so that they are not in the plot
    mask = res_e == 0
    res_e[mask] = np.nan

    return bins, res_e


def energy_bias(
        true_energy, reco_energy, bins=None, 
        x_axis_true_energy=True):
    """
    Evaluates energy bias per bin in energy

    Parameters
    ----------
    true_energy: astropy.units.quantity.Quantity
    reco_energy: astropy.units.quantity.Quantity
    bins: astropy.units.quantity.Quantity
        Energy bin edges
    x_axis_true_energy: bool
        If True energy axis is binned in true
        energy. Reconstructed energy is used otherwise

    Returns
    -------
    bins: astropy.units.quantity.Quantity
    bias_e: numpy.ndarray

    """

    bias_e = []

    for i in range(len(bins)-1):
        if x_axis_true_energy:
            mask = (true_energy > bins[i]) & (true_energy < bins[i + 1])
        else:
            mask = (reco_energy > bins[i]) & (reco_energy < bins[i + 1])
        bias_e.append(ctaplot.ana.relative_bias(true_energy[mask], reco_energy[mask], relative_scaling_method='s1'))

    bias_e = np.array(bias_e, dtype=float)
    # replace zeros with Nans so that they are not in the plot
    mask = bias_e == 0
    bias_e[mask] = np.nan

    return bins, np.array(bias_e)


def angular_resolution_per_energy(
        true_alt, reco_alt, true_az, reco_az, energy,
        percentile=68.27, confidence_level=0.95, 
        bias_correction=False, bins=None):
    """
    Evaluates angular resolution per bin in energy

    Parameters
    ----------
    true_alt: astropy.units.quantity.Quantity
    reco_alt: astropy.units.quantity.Quantity
    true_az: astropy.units.quantity.Quantity
    reco_az: astropy.units.quantity.Quantity
    energy: astropy.units.quantity.Quantity
    percentile: float
    confidence_level: float
    bias_correction: bool
    bins: astropy.units.quantity.Quantity
        Energy bin edges

    Returns
    -------
    bins: astropy.units.quantity.Quantity
    res_q: astropy.units.quantity.Quantity

    """

    if not len(reco_alt) == len(reco_az) == len(energy) > 0:
        raise ValueError("reco_alt, reco_az and true_energy must have the same length")

    res = []

    logging.info(f"Out of {len(reco_alt)}, in {sum(~(np.isfinite(reco_alt) & np.isfinite(reco_az)))} events the altitude or azimuth was not reconstructed")

    for i, e in enumerate(bins[:-1]):
        mask = (energy > bins[i]) & (energy <= bins[i + 1]) & np.isfinite(reco_alt) & np.isfinite(reco_az)
        res.append(ctaplot.ana.angular_resolution(true_alt[mask], reco_alt[mask], true_az[mask], reco_az[mask],
                                            percentile=percentile,
                                            confidence_level=confidence_level,
                                            bias_correction=bias_correction,
                                            )
                   )

    res_q = u.Quantity(res)

    # replace zeros with Nans so that they are not in the plot
    mask = res_q == 0 * u.rad
    res_q[mask] = np.nan

    return bins, res_q.to(u.deg)


def energy_resolution(
        dl2_gamma, e_bins=None, outdir=None, 
        telescope=None, save_fig=False, save_hdf=False, 
        x_axis_true_energy=True):
    """
    Evaluates energy resolution and bias.

    Parameters
    ----------
    dl2_gamma: pandas.DataFrame
        DL2 MC table for gammas (after some optional 
        gammaness cutting). True and reconstructed energy
        is expected to be in TeV
    e_bins: astropy.units.quantity.Quantity
        Energy bin edges
    outdir: string
    telescope: string
    save_fig: bool
        If True some figures are stored
    save_hdf: bool
        If True hdf files with performance
        tables are stored
    x_axis_true_energy: bool
        If True energy axis is binned in true
        energy. Reconstructed energy is used otherwise

    Returns
    -------

    """

    # energy resolution
    # We cannot use ctaplot.plot_energy_resolution(), because in this function I am not able to 
    # set custom energy bins, it's somehow messed up
    # NOTE: Be carefull with the units! - must be revisited for proper unit conversion!
    e_bin, e_res = energy_resolution_per_energy(
        dl2_gamma.true_energy.values * u.TeV, 
        dl2_gamma.reco_energy.values * u.TeV, 
        bins=e_bins, 
        x_axis_true_energy=True
        )

    # Energy bias = np.median((reco - true) / np.abs(true)
    e_bin, e_bias = energy_bias(
        dl2_gamma.true_energy.values * u.TeV, 
        dl2_gamma.reco_energy.values * u.TeV, 
        bins=e_bins,
        x_axis_true_energy=True
        )

    if save_fig:

        logging.info('Energy resolution: saving images..')
        energy_center = ctaplot.ana.logbin_mean(e_bins)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))

        axes[0, 0].set_ylabel(r"$(\Delta energy/energy)_{68}$")
        if x_axis_true_energy:
            axes[0, 0].set_xlabel(rf'$E_T$ [{energy_center.unit.to_string("latex")}]')
        else:
            axes[0, 0].set_xlabel(rf'$E_R$ [{energy_center.unit.to_string("latex")}]')
        axes[0, 0].set_xscale('log')
        axes[0, 0].set_title('Energy resolution')

        axes[0, 0].errorbar(energy_center, e_res[:, 0],
                    xerr=(energy_center - e_bin[:-1], e_bin[1:] - energy_center),
                    yerr=(e_res[:, 0] - e_res[:, 1], e_res[:, 2] - e_res[:, 0]),
                    fmt='o'
                    )

        axes[0, 0].set_xlim([min(e_bins.value), max(e_bins.value)])
        axes[0, 0].set_ylim([0, 0.7])
        axes[0, 0].grid(True, which='both')

        axes[1, 0].set_ylabel(r"bias (median($E_{reco}/E_{true}$ - 1)")
        if x_axis_true_energy:
            axes[1, 0].set_xlabel(rf'$E_T$ [{energy_center.unit.to_string("latex")}]')
        else:
            axes[1, 0].set_xlabel(rf'$E_R$ [{energy_center.unit.to_string("latex")}]')
        axes[1, 0].set_xscale('log')
        axes[1, 0].set_title('Energy bias')

        axes[1, 0].errorbar(energy_center, e_bias, xerr=(energy_center - e_bin[:-1], e_bin[1:] - energy_center), fmt='o')
        axes[1, 0].grid(True, which='both')

        ctaplot.plot_migration_matrix(dl2_gamma.true_energy.apply(np.log10),
                                    dl2_gamma.reco_energy.apply(np.log10),
                                    ax=axes[0, 1],
                                    colorbar=True,
                                    xy_line=True,
                                    hist2d_args=dict(norm=mpl.colors.LogNorm(), range=np.array([(-1, 2.5), (-1,2.5)])),
                                    line_args=dict(color='black'),
                                    )

        axes[0, 1].set_xlabel('log(mc energy/[TeV])')
        axes[0, 1].set_ylabel('log(reco energy/[TeV])')
        axes[0, 0].set_title("")
        axes[0, 0].label_outer()
        axes[1, 0].set_title("")
        axes[1, 0].set_ylabel("Energy bias")
        axes[1, 0].set_xlim([min(e_bins.value), max(e_bins.value)])
        axes[1, 0].set_ylim([-0.3, 0.3])
        axes[1, 1].remove()

        fig.tight_layout()
        fig.savefig(outdir + '/energy_resolution_'+telescope+'.png', dpi=200)

    if save_hdf:

        logging.info(f'Energy resolution: saving tables for {telescope}..')

        e_bins_t = QTable(data=e_bins[..., np.newaxis], names=['energy_bins'])

        data = e_res
        names = ['energy_res', 'energy_res_err_lo', 'energy_res_err_hi']
        data = np.append(data, e_bias[..., np.newaxis], axis=1)
        names.append('energy_bias')
        res_t = Table(data=data, names=names)

        outfile = outdir + '/energy_resolution_'+telescope+'.h5'
        write_table_hdf5(e_bins_t, outfile, path='bins', overwrite=True, append=True, serialize_meta=True)
        write_table_hdf5(res_t, outfile, path='res', overwrite=True, append=True)


def angular_resolution(
        dl2_gamma, e_bins=None, outdir=None, 
        telescope=None, save_fig=False, save_hdf=False, 
        x_axis_true_energy=True, axes_sky=True):
    """
    Evaluates angular resolution.

    Parameters
    ----------
    dl2_gamma: pandas.DataFrame
        DL2 MC table for gammas (after some optional 
        gammaness cutting)
    e_bins: astropy.units.quantity.Quantity
        Energy bin edges
    outdir: string
    telescope: string
    save_fig: bool
        If True some figures are stored
    save_hdf: bool
        If True hdf files with performance
        tables are stored
    x_axis_true_energy: bool
        If True energy axis is binned in true
        energy. Reconstructed energy is used otherwise
    axes_sky: bool
        Only affets plotting. If True, axes of 
        the 2D histogram are in alt/az, or in
        x/y (camera frame) otherwise.

    Returns
    -------

    """

    # Angular resolution
    # NOTE: Be carefull with the units! - must be revisited for proper unit conversion!
    # True or reco energy? - in lstchain, there is reco_energy, ctaplot is ambivalent
    if x_axis_true_energy:
        energy_x = dl2_gamma.true_energy.values *  u.TeV
    else:
        energy_x = dl2_gamma.reco_energy.values *  u.TeV
    e_bin, ang_res = angular_resolution_per_energy(dl2_gamma.true_alt.values * np.pi / 180. * u.rad,
                                                dl2_gamma.reco_alt.values * np.pi / 180. *  u.rad,
                                                dl2_gamma.true_az.values * np.pi / 180. *  u.rad,
                                                dl2_gamma.reco_az.values * np.pi / 180. *  u.rad,
                                                energy_x,
                                                bins=e_bins)

    if save_fig:

        logging.info('Angular resolution: saving images..')
        energy_center = ctaplot.ana.logbin_mean(e_bins)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))

        # Reconstructed point gamma positions
        if axes_sky:
            hist1 = axes[0, 0].hist2d(dl2_gamma['reco_az'], dl2_gamma['reco_alt'], bins=50, norm=mpl.colors.LogNorm())
            axes[0, 0].set_xlabel('az [deg]')
            axes[0, 0].set_ylabel('alt [deg]')

            dl2_gamma_azimuth_transformed = (dl2_gamma['true_az'] +  180) % 360 - 180
            if abs(dl2_gamma['true_az_tel'].iloc[0]) <= 90:
                logging.info('Using transformed azimuth to plot angular resolution.')
                axes[0, 0].set_xlim(dl2_gamma_azimuth_transformed.iloc[0]-5/np.cos(np.deg2rad(dl2_gamma['true_alt'].iloc[0])), dl2_gamma_azimuth_transformed.iloc[0]+5/np.cos(np.deg2rad(dl2_gamma['true_alt'].iloc[0])))
            else:
                logging.info('Using non-transformed azimuth to plot angular resolution.')
                axes[0, 0].set_xlim(dl2_gamma['true_az'].iloc[0]-5/np.cos(np.deg2rad(dl2_gamma['true_alt'].iloc[0])), dl2_gamma['true_az'].iloc[0]+5/np.cos(np.deg2rad(dl2_gamma['true_alt'].iloc[0])))

            axes[0, 0].set_ylim(dl2_gamma['true_alt'].iloc[0]-5, dl2_gamma['true_alt'].iloc[0]+5)
        elif not axes_sky and 'reco_src_x' in dl2_gamma.keys():
            hist1 = axes[0, 0].hist2d(dl2_gamma['reco_src_x'], dl2_gamma['reco_src_y'], bins=50, range=np.array([(-1, 1), (-1, 1)]), norm=mpl.colors.LogNorm())
            axes[0, 0].set_xlabel('x [m]')
            axes[0, 0].set_ylabel('y [m]')
        else:
            raise Exception('ERROR: cannot plot reconstructed image in camera frame, because reco_src_x/y are not in DL2 file.')

        cbar = plt.colorbar(hist1[3], ax=axes[0, 0], orientation='vertical')
        cbar.set_label('N of events')
        axes[0, 0].grid()

        # Theta2 plot
        # NOTE: Be carefull with the units! - must be revisited for proper unit conversion!
        ctaplot.plot_theta2(dl2_gamma.true_alt.values * u.deg,
                    dl2_gamma.reco_alt.values * u.deg,
                    dl2_gamma.true_az.values * u.deg,
                    dl2_gamma.reco_az.values * u.deg,
                    ax=axes[0, 1],
                    bins=100,
                    range=(0, 0.5),
                    )
        axes[0, 1].grid(True, which='both')

        # Angular resolution is traditionally presented in degrees
        ang_res = ang_res.to(u.deg)

        axes[1, 0].set_ylabel('Angular Resolution [deg]')
        if x_axis_true_energy:
            axes[1, 0].set_xlabel(rf'$E_T$ [{energy_center.unit.to_string("latex")}]')
        else:
            axes[1, 0].set_xlabel(rf'$E_R$ [{energy_center.unit.to_string("latex")}]')
        axes[1, 0].set_xscale('log')
        axes[1, 0].set_title('Angular resolution')
        axes[1, 0].grid(True, which='both')
        axes[1, 0].errorbar(energy_center, 
                        ang_res[:, 0], 
                        xerr=(energy_center - e_bin[:-1], e_bin[1:] - energy_center),
                        yerr=(ang_res[:, 0] - ang_res[:, 1], ang_res[:, 2] - ang_res[:, 0]),
                        fmt='o')
        axes[1, 0].set_xlim([min(e_bins.value), max(e_bins.value)])
        axes[1, 0].set_ylim([0, 1])

        axes[1, 1].remove()

        fig.tight_layout()
        fig.savefig(outdir + '/angular_resolution_'+telescope+'.png', dpi=200)

    if save_hdf:

        logging.info(f'Angular resolution: saving tables for {telescope}..')

        e_bins_t = QTable(data=e_bins[..., np.newaxis], names=['energy_bins'])

        data = ang_res
        names = ['angular_res', 'angular_res_err_lo', 'angular_res_err_hi']
        res_t = Table(data=data, names=names)

        outfile = outdir + '/angular_resolution_'+telescope+'.h5'
        write_table_hdf5(e_bins_t, outfile, path='bins', overwrite=True, append=True, serialize_meta=True)
        write_table_hdf5(res_t, outfile, path='res', overwrite=True, append=True)


def roc_curve(
        gh_testing_dataset, e_bins=None, 
        outdir=None, telescope=None, save_fig=False, 
        save_hdf=False):
    """
    Evaluates ROC for gamma/hadron separation
    performance.

    Parameters
    ----------
    gh_testing_dataset: pandas.DataFrame
        DL2 MC mixed table of protons and gammas
    e_bins: astropy.units.quantity.Quantity
        Energy bin edges
    outdir: string
    telescope: string
    save_fig: bool
        If True some figures are stored
    save_hdf: bool
        If True hdf files with performance
        tables are stored

    Returns
    -------

    """

    if save_fig:

        logging.info('ROC curve: saving images..')

        fig, ax = plt.subplots(1, 2, figsize=(12, 6))

        ctaplot.plot_roc_curve_gammaness(gh_testing_dataset.true_shower_primary_id, gh_testing_dataset.gammaness,
                                              ax=ax[0],
                                              )
        ax[0].set_title("gamma ROC curve")
        ax[0].set_xlabel("gamma false positive rate")
        ax[0].set_ylabel("gamma true positive rate")
        ax[0].grid(True, which='both')
        ax[0].set_xlim([0, 1])
        ax[0].set_ylim([0, 1])
    
        # gammaness distribution of protons and gammas
        mask_gamma = gh_testing_dataset.true_shower_primary_id == 0
        mask_proton = gh_testing_dataset.true_shower_primary_id == 101

        ax[1].hist(gh_testing_dataset[mask_gamma].gammaness, range=(0, 1), bins = 100, density=True, label='gammas', alpha=0.5)
        ax[1].hist(gh_testing_dataset[mask_proton].gammaness, range=(0, 1), bins = 100, density=True, label='protons', alpha=0.5)
        ax[1].set_xlabel("gammaness")
        ax[1].legend()
        ax[1].set_xlim([0, 1])
        ax[1].grid(True, which='both')

        fig.tight_layout()
        fig.savefig(outdir + '/roc_curve_'+telescope+'.png', dpi=200)

        # ROC and gammaness in reco energy bins
        energy = gh_testing_dataset['reco_energy']
        for i, e in enumerate(e_bins[:-1]):

            mask_e = (energy > e_bins[i]) & (energy <= e_bins[i + 1])
            if (sum(mask_e) > 10) & (sum(gh_testing_dataset[mask_e].true_shower_primary_id == 0) > 0) & (sum(gh_testing_dataset[mask_e].true_shower_primary_id == 101) > 0):

                emin = round(e_bins[i].value, 2)
                emax = round(e_bins[i+1].value, 2)

                fig, ax = plt.subplots(1, 2, figsize=(12, 6))

                ctaplot.plot_roc_curve_gammaness(gh_testing_dataset[mask_e].true_shower_primary_id, gh_testing_dataset[mask_e].gammaness,
                                                    ax=ax[0],
                                                    )
                #ax[0].set_title("gamma ROC curve, E = [%.2f, %.2f]", emin, emax)
                ax[0].set_xlabel("gamma false positive rate")
                ax[0].set_ylabel("gamma true positive rate")
                ax[0].grid(True, which='both')
                ax[0].set_xlim([0, 1])
                ax[0].set_ylim([0, 1])

                ax[1].hist(gh_testing_dataset[mask_e & mask_gamma].gammaness, range=(0, 1), bins = 100, density=True, label='gammas', alpha=0.5)
                ax[1].hist(gh_testing_dataset[mask_e & mask_proton].gammaness, range=(0, 1), bins = 100, density=True, label='protons', alpha=0.5)
                ax[1].set_xlabel("gammaness")
                ax[1].legend()
                ax[1].set_xlim([0, 1])
                ax[1].grid(True, which='both')

                outdir_roc = outdir + '/roc_bins/'
                check_outdir(outdir_roc)
                fig.tight_layout()
                fig.savefig(outdir_roc + '/roc_curve_'+telescope+'_'+str(emin)+'_'+str(emax) + '.png', dpi=200)
                plt.close()

    if save_hdf:

        logging.info(f'ROC curve: saving tables for {telescope}..')

        fpr, tpr, thresholds = metrics.roc_curve(gh_testing_dataset.true_shower_primary_id,
                                             gh_testing_dataset.gammaness,
                                             pos_label=0,
                                             sample_weight=None,
                                             drop_intermediate=True,
                                             )

        data = np.array(np.column_stack((fpr, tpr, thresholds)))

        names = ['false_positive_rate', 'true_positive_rate', 'threshold']
        roc = Table(data=data, names=names)
        outfile = outdir + '/roc_curve_'+telescope+'.h5'
        write_table_hdf5(roc, outfile, path='roc', overwrite=True, append=True, serialize_meta=True)


class irf_maker:

    def __init__(self,
                 config_filename    = None,
                 mc_gamma_filename  = None,
                 mc_proton_filename = None,
                 mc_tel_setup       = None,
                 point_like_offset  = None,
                 output_dir = './data/',
                 gammaness_cuts = None,
                 true_energy_scaling = False,
                 event_class = None,
                 ):

        # This is arbitrary as the background is normalized in the analysis
        self.bg_obstime = 50*u.h

        self.config_filename    = config_filename
        self.mc_gamma_filename  = mc_gamma_filename
        self.mc_proton_filename = mc_proton_filename
        self.output_dir = output_dir

        self.config = load_config(self.config_filename)

        # True energy scaling
        if true_energy_scaling:
            self.scaling_factor = float(self.config["analysis"]["true_energy_scaling_factor"])
            logging.warning('True energies in IRFs scaled by a factor of %f.', self.scaling_factor)
        else: self.scaling_factor = 1.

        if gammaness_cuts is None:
            self.gammaness_cut = self.config['analysis']['global_gammaness_cut']
            self.gammaness_cut_tag = self.gammaness_cut
        else:
            self.gammaness_cut_tag = 'energydep'
            self.gammaness_cut = read_table_hdf5(gammaness_cuts, path='gammaness_cuts')

        self.point_like_offset = point_like_offset
        v1 = sst1mpipe.__version__.split('.')[0]
        v2 = sst1mpipe.__version__.split('.')[1]
        self.pipeline_version = '{}_{}'.format(v1,v2)

        self.hdu_list = [fits.PrimaryHDU()]

        if mc_tel_setup[-1]=='1':
            tel_setup='tel_021'
        elif mc_tel_setup[-1]=='2':
            tel_setup='tel_022'
        else:
            tel_setup=mc_tel_setup
        self.tel_setup = tel_setup

        logging.info("Making IRFs for telescope {}".format(self.tel_setup))

    ##################################### "dl2 tables" #####################################

        sim_info = get_mc_info(mc_gamma_filename)
        self.sim_info = SimulatedEventsInfo(
                               n_showers      = sim_info.n_showers, 
                               energy_min     = sim_info.energy_min, 
                               energy_max     = sim_info.energy_max, 
                               max_impact     = sim_info.max_impact, 
                               spectral_index = sim_info.spectral_index, 
                               viewcone       = sim_info.viewcone) ## min & max in latter version

        dl2_mc_gamma  = load_dl2_sst1m(mc_gamma_filename,
                                       config=self.config,
                                       tel=mc_tel_setup)

        dl2_mc_proton = load_dl2_sst1m(mc_proton_filename,
                                       config=self.config,
                                       tel=mc_tel_setup)
        sim_info_proton = get_mc_info(mc_proton_filename)

        if true_energy_scaling:
            dl2_mc_gamma['true_energy'] *= self.scaling_factor
            dl2_mc_proton['true_energy'] *= self.scaling_factor
            self.sim_info.energy_min *= self.scaling_factor
            self.sim_info.energy_max *= self.scaling_factor
            sim_info_proton.energy_min *= self.scaling_factor
            sim_info_proton.energy_max *= self.scaling_factor

        dl2_mc_proton = get_weights(dl2_mc_proton, 
                                    mc_info  = sim_info_proton, 
                                    obs_time = self.bg_obstime, 
                                    target_spectrum = DAMPE_P_He_SPECTRUM)
        
        if event_class is not None:
            logging.info("Making IRFs for event class {}".format(event_class))
            dl2_mc_gamma_c = dl2_mc_gamma[dl2_mc_gamma['event_type'] == event_class]
            dl2_mc_proton_c = dl2_mc_proton[dl2_mc_proton['event_type'] == event_class]
            logging.info("{} gammas of event class {}".format(len(dl2_mc_gamma_c), event_class))
            logging.info("{} protons of event class {}".format(len(dl2_mc_proton_c), event_class))
        else:
            dl2_mc_gamma_c = dl2_mc_gamma
            dl2_mc_proton_c = dl2_mc_proton

        self.gamma_event_dict  = self.table_to_selectedEvt_dict(dl2_mc_gamma_c)
        self.proton_event_dict = self.table_to_selectedEvt_dict(dl2_mc_proton_c)

        logging.info("{} gamma events after selection".format(len(self.gamma_event_dict['true_energy'])))
        logging.info("{} proton events after selection".format(len(self.proton_event_dict['true_energy'])))

        ######################## Set outdir and out filenames ########################
        tel_ze, tel_az = get_avg_pointing(dl2_mc_gamma_c)
        self.zenith_angle      = int(round(tel_ze/10.)*10)
        self.azimuth      = int(round(tel_az/10.)*10)

        # /data/... is mandatory scheme. if not used, the hdu indexer cannot merge photon lists with IRFs
        # TODO: deal with the NSB bin if needed, removing NSB from the IRF name for now
        #self.outdir = output_dir + '/data/sst1m_{}/{}/bcf/ze{}_az{}_nsb100_gc{}/'.format(self.tel_setup,
        #                                                            self.pipeline_version,
        #                                                            self.zenith_angle,
        #                                                            self.azimuth,
        #                                                            self.gammaness_cut_tag)
        if event_class is not None:
            ec_str = "_ec{}".format(event_class)
        else: ec_str = ""
        self.outdir = output_dir + '/data/sst1m_{}/{}/bcf/ze{}_az{}_gc{}{}/'.format(self.tel_setup,
                                                                    self.pipeline_version,
                                                                    self.zenith_angle,
                                                                    self.azimuth,
                                                                    self.gammaness_cut_tag,
                                                                    ec_str)

        if point_like_offset is not None :
            ptlk_str = "_pointlike_{}deg".format(point_like_offset)
        else:
            ptlk_str = ""
        self.out_fits_filename = "SST1M_{}_Zen{}deg_gcut{}{}{}_irfs.fits".format(self.tel_setup,
                                                                               self.zenith_angle,
                                                                               self.gammaness_cut_tag,
                                                                               ptlk_str,
                                                                               ec_str)
        
    ##################################################################################
    ##################################### "BINS" #####################################
        
        emin = self.config['analysis']["log_energy_min_tev"]
        emax = self.config['analysis']["log_energy_max_tev"]
        nbinse = self.config['analysis']["n_energy_bins"]
        self.reco_energy_bins = np.logspace(emin, emax, nbinse) * u.TeV
        self.true_energy_bins = np.geomspace(
            sim_info.energy_min.to_value(u.TeV), 
            sim_info.energy_max.to_value(u.TeV), 
            nbinse
            ) * u.TeV
        
        nbins_migration = self.config['analysis']["n_bins_migration"]
        self.e_bins_edisp = np.geomspace(
            sim_info.energy_min.to_value(u.TeV), 
            sim_info.energy_max.to_value(u.TeV), 
            nbins_migration
            ) * u.TeV

        # From pyirf docu: Bin edges in relative deviation (not bins in energy!)
        self.energy_migration_bins = np.linspace(0.1, 4, nbins_migration)

        fov_offset_max = self.config['analysis']["fov_offset_max_deg"]
        fov_offset_bg_max = self.config['analysis']["fov_offset_bg_max_deg"]
        source_offset_max = self.config['analysis']["source_offset_max_deg"]
        nbins_fov = self.config['analysis']["nbins_fov_offset"]
        nbins_fov_bg = self.config['analysis']["nbins_fov_offset_bg"]
        nbins_offset = self.config['analysis']["nbins_source_offset"]

        self.fov_offset_bins_bg    = np.linspace(0, fov_offset_bg_max, nbins_fov_bg) * u.deg
        self.source_offset_bins =    np.linspace(0, source_offset_max, nbins_fov_bg)* u.deg

        if point_like_offset is not None:
            self.fov_offset_bins    = np.linspace(max(0,point_like_offset-1),point_like_offset+1, 2) * u.deg
            self.point_like = True
        else:
            self.fov_offset_bins    = (np.linspace(0, fov_offset_max, nbins_fov) * u.deg)
            self.point_like = False

################# table_to_selectedEvt_dict ###########

    def table_to_selectedEvt_dict(self, dl2_data):

        if isfloat(self.gammaness_cut):
            logging.info('Global gammaness cut {} applied.'.format(self.gammaness_cut))
            mask_gg = dl2_data['gammaness'] > self.gammaness_cut
        else:
            logging.info('Energy dependent gammaness cut applied.')

            mask_gg = evaluate_binned_cut(
                dl2_data["gammaness"],
                dl2_data["reco_energy"] * u.TeV,
                self.gammaness_cut,
                operator.ge,
            )
        # JJ: .copy() complains: ValueError: values whose keys begin with an uppercase char must be Config instances: 'MC_correction_for_PDE', True
        # I changed the config param to lower case, let's see what happens in the next versions
        dl2_selected = dl2_data[mask_gg] #.copy()

        # event selection is performed authomaticaly, if you provide load_dl2_sst1m with a config file
        # dl2_selected = event_selection(dl2_selected, config=self.config)
        tel_altaz = AltAz( alt=dl2_selected["true_alt_tel"]*u.deg, az=dl2_selected["true_az_tel"]*u.deg)
        evt_true_altaz = AltAz(alt=dl2_selected["true_alt"]*u.deg, az=dl2_selected["true_az"]*u.deg)
        evt_reco_altaz = AltAz(alt=dl2_selected["reco_alt"]*u.deg, az=dl2_selected["reco_az"]*u.deg)

        event_dict = {}
        for key in ['true_energy','reco_energy']:
            event_dict[key] = np.array(dl2_selected[key]) * u.TeV
        event_dict['reco_source_fov_offset'] = tel_altaz.separation(evt_reco_altaz).to(u.deg)
        event_dict['true_source_fov_offset'] = tel_altaz.separation(evt_true_altaz).to(u.deg)
        event_dict['theta'] = evt_reco_altaz.separation(evt_true_altaz).to(u.deg)
        try:
            event_dict['weight'] = dl2_selected['weight']
        except:
            pass
        return event_dict
###################################################
###################################################

    def check_irf_exist(self):
        pass
    
    ############## AEFF ##############
    def make_aeff_irf(self):
                
        if self.point_like :
            
            logging.info("Making Effective area (point-like).")
            aeff = effective_area_per_energy(selected_events  = Table(self.gamma_event_dict),
                                             simulation_info  = self.sim_info,
                                             true_energy_bins = self.true_energy_bins)
            hdu_aeff=create_aeff2d_hdu(effective_area   = aeff[:, np.newaxis], 
                                       true_energy_bins = self.true_energy_bins,
                                       fov_offset_bins  = self.fov_offset_bins,
                                       point_like       = self.point_like)

        else :
            logging.info("Making Effective area (full enclosure).")
            aeff = effective_area_per_energy_and_fov(selected_events  = self.gamma_event_dict,
                                                     simulation_info  = self.sim_info,
                                                     true_energy_bins = self.true_energy_bins,
                                                     fov_offset_bins  = self.fov_offset_bins)

            hdu_aeff=create_aeff2d_hdu(effective_area   = aeff, 
                                       true_energy_bins = self.true_energy_bins,
                                       fov_offset_bins  = self.fov_offset_bins,
                                       extname          = "EFFECTIVE AREA",
                                       point_like       = self.point_like)
        return hdu_aeff

    ############## EDISP ##############
    def make_edisp_irf(self):
        
        logging.info("Making Energy dispersion matrix.")
        edisp = energy_dispersion(selected_events  =self.gamma_event_dict, 
                                  true_energy_bins =self.e_bins_edisp, 
                                  fov_offset_bins  =self.fov_offset_bins, 
                                  migration_bins   =self.energy_migration_bins)

        hdu_edisp = create_energy_dispersion_hdu(energy_dispersion = edisp,
                                                 true_energy_bins  = self.e_bins_edisp,
                                                 migration_bins    = self.energy_migration_bins,
                                                 fov_offset_bins   = self.fov_offset_bins,
                                                 extname           = "ENERGY DISPERSION",
                                                 point_like        = self.point_like)
        return hdu_edisp


    ############## PSF ##############
    def make_psf_irf(self):

        logging.info("Making Point Spread Function.")
        psf = psf_table(events             = self.gamma_event_dict, 
                        true_energy_bins   = self.true_energy_bins, 
                        source_offset_bins = self.source_offset_bins, 
                        fov_offset_bins    = self.fov_offset_bins
                        )
        hdu_psf=create_psf_table_hdu(psf                = psf,
                                     true_energy_bins   = self.true_energy_bins,
                                     source_offset_bins = self.source_offset_bins,
                                     fov_offset_bins    = self.fov_offset_bins,
                                     extname            = "POINT SPREAD FUNCTION",
                                     point_like         = self.point_like)
        return hdu_psf


    ############## BKG ##############
    def make_bkg_irf(self):
        
        logging.info("Making Background model.")
        bg_2d = background_2d(events           = self.proton_event_dict,
                              reco_energy_bins = self.reco_energy_bins,
                              fov_offset_bins  = self.fov_offset_bins_bg,
                              t_obs            = self.bg_obstime)

        hdu_bg2d = create_background_2d_hdu(background_2d    = bg_2d,
                                            reco_energy_bins = self.reco_energy_bins,
                                            fov_offset_bins  = self.fov_offset_bins_bg,
                                            extname          = "BACKGROUND",
                                            point_like       = self.point_like)
        #hdu_bg2d.header['HDUCLAS4']='3D'
        return hdu_bg2d

    ###########################


    def make_all_irfs(self, save = True):

        self.hdu_list.append(self.make_aeff_irf())
        self.hdu_list.append(self.make_psf_irf())
        self.hdu_list.append(self.make_edisp_irf())
        self.hdu_list.append(self.make_bkg_irf())

        if save :
            check_outdir(self.outdir)
            fits.HDUList(self.hdu_list).writeto(self.outdir+self.out_fits_filename,
                                                 overwrite=True)
            # Config file cannot be placed in the directory with IRFs, because otherwise HDU indexer 
            # is not able to do the match
            shutil.copy(self.config_filename,self.output_dir)
            logging.info("IRFs stored in: {}".format(self.outdir))

def isfloat(num):
    try:
        float(num)
        return True
    except TypeError:
        return False
