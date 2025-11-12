import matplotlib.pyplot as plt
import numpy as np
from gammapy.stats import WStatCountsStatistic
import astropy.units as u
import itertools 
from astropy.convolution import convolve
from astropy.convolution.kernels import Gaussian2DKernel
from astropy.visualization.wcsaxes import SphericalCircle
import ctaplot
from sst1mpipe.performance.spectra import *


def plot_count_maps(
        data, bins_raw=20, bins_conv=100, 
        conv_kernel_deg=0.2, range_deg=4, 
        source=None, theta2_cut=None, 
        wobbles=[], off_radec_w=[]):
    """
    Plots count maps (not bkg subtracted). Raw 2D historam and
    the same but convoluted with a gaussian kernel of given
    size.

    Parameters
    ----------
    data: pandas.DataFrame with DL3 table
    bins_raw: int
        Bins in raw 2D histogram
    bins_conv: int
        Bins in convaluted 2D histogram
    conv_kernel_deg: float
        Size of the convolution kernel
    range_deg: float
    source: astropy.coordinates.SkyCoord
    theta2_cut: float
    wobbles: list of astropy.coordinates.SkyCoord
    off_radec_w: list of astropy.coordinates.SkyCoord
        Ra/dec coordinates of the OFF regions
    Returns
    -------

    """

    range_d = [
                       [source.ra.value-range_deg/2/np.cos(np.deg2rad(source.dec.value)), 
                        source.ra.value+range_deg/2/np.cos(np.deg2rad(source.dec.value))], 
                       [source.dec.value-range_deg/2, 
                        source.dec.value+range_deg/2]
                   ]
    
    fig, ax = plt.subplots(1, 2, figsize=(17,7))

    h = ax[0].hist2d(data.reco_ra, 
                   data.reco_dec, 
                   bins=bins_raw, 
                   range=range_d
                  )

    ax[0].set_xlabel('RA [deg]')
    ax[0].set_ylabel('DEC [deg]')

    # Convolution with a gaussian kernel
    image_range_deg = np.array(range_d)

    n_pix_deg = bins_conv / range_deg  # px / deg
    sigma_px = conv_kernel_deg * n_pix_deg    # px
    kernel_x_size = int(5*sigma_px)

    if kernel_x_size % 2 == 0:  # kernel size must be odd
         kernel_x_size =  kernel_x_size + 1
    kernel_y_size = kernel_x_size

    kernel = Gaussian2DKernel(sigma_px, mode='oversample', factor=100, x_size=kernel_x_size, y_size=kernel_y_size).array
    heatmap, xedges, yedges = np.histogram2d(
                            data.reco_ra, 
                            data.reco_dec, 
                            bins=bins_conv, 
                            range=image_range_deg
                            )
    image = convolve(heatmap, kernel)

    # Convoluted image
    edge = int(0.05*bins_conv) # We do not plot the borders affected by convolution artefacts
    if edge > 0:
        img = ax[1].pcolor(xedges[edge:-edge], yedges[edge:-edge], image[edge:-edge, edge:-edge].T)
    else:
        img = ax[1].pcolor(xedges, yedges, image.T)

    ax[1].set_xlabel('RA [deg]')
    ax[1].set_ylabel('DEC [deg]')
    ax[1].grid(alpha=0.7)
    
    #fig.colorbar(h[3], ax=ax[0])
    #fig.colorbar(img, ax=ax[1])

    if source is not None:
        on1 = SphericalCircle((source.ra, source.dec), theta2_cut**0.5 * u.deg, resolution=100, fill=False, color='red', linewidth=2, alpha=0.7)
        ax[0].add_patch(on1)
        # If not used like that: "ValueError: Can not reset the axes.  You are probably trying to re-use an artist in more than one Axes which is not supported"
        on2 = SphericalCircle((source.ra, source.dec), theta2_cut**0.5 * u.deg, resolution=100, fill=False, color='red', linewidth=2, alpha=0.7)
        ax[1].add_patch(on2)

    if len(wobbles) > 0:
        for wobble in wobbles:
            ax[0].plot(wobble.ra.value, wobble.dec.value, 'x', color='red', markersize=10)
            ax[0].text(wobble.ra.value, wobble.dec.value + range_deg/50, wobble.info.name, color='red', fontsize=15)
            ax[1].plot(wobble.ra.value, wobble.dec.value, 'x', color='red', markersize=10)
            ax[1].text(wobble.ra.value, wobble.dec.value + range_deg/50, wobble.info.name, color='red', fontsize=15)

    if len(off_radec_w) > 0:
        for off_w in off_radec_w:
            off1 = SphericalCircle((off_w.ra, off_w.dec), theta2_cut**0.5 * u.deg, resolution=100, fill=False, color='white', linewidth=3)
            ax[0].add_patch(off1)
            # If not used like that: "ValueError: Can not reset the axes.  You are probably trying to re-use an artist in more than one Axes which is not supported"
            off2 = SphericalCircle((off_w.ra, off_w.dec), theta2_cut**0.5 * u.deg, resolution=100, fill=False, color='white', linewidth=3)
            ax[1].add_patch(off2)


def cycle(iterable):
    for item in itertools.cycle(iterable):
        yield item


def plot_sigma_time(data, sigma, times, obsid, alphas, nights):
    """
    Plots sigma and alpha evolution with time. Alpha here stands for
    bkg normalization factor. It provides also a plot with per-night
    color of sigma=f(time) for easy identification of low-quality datasets.

    Parameters
    ----------
    data: astropy.table.table.Table
    sigma: list
    times: list
    obsid: list
    alphas: list
    nights: list of ints

    Returns
    -------

    """

    fig, ax = plt.subplots(1, 3, figsize=(25,7))
    ax[0].plot(np.sqrt(times), sigma, 'b.')
    ax[0].plot([0, np.sqrt(times[-1])], [0, max(sigma)], 'y-')
    ax[0].grid()
    ax[0].set_xlabel(r'$\sqrt{t} \, [\sqrt{s}]$')
    ax[0].set_ylabel('$\sigma$')

    ax[1].plot(times, alphas, 'b.')
    ax[1].grid()
    ax[1].set_xlabel(r'$t \, [s]$')
    ax[1].set_ylabel(r'$\alpha$')
    
    times = np.array(times)
    sigma = np.array(sigma)
    obsid = np.array(obsid)

    colors = ['r', 'b', 'g', 'y', 'm']

    for i, c in zip(range(len(nights)), cycle(colors)):
        if i < len(nights)-1:
            mask = (obsid >= nights[i]) & (obsid <= nights[i+1])
        else:
            mask = (obsid >= nights[i])
        ax[2].plot(times[mask]/3600, sigma[mask], '.', color=c)
        
        ax[2].text(times[mask][0]/3600, sigma[mask][0]-1, str(nights[i])[4:8], fontsize=12, color=c, fontweight='bold')

    ax[2].plot(times/3600, (max(sigma) / np.sqrt(times[-1]/3600)) * np.sqrt(times/3600), 'y-')
    ax[2].grid()
    ax[2].set_xlabel('t [h]')
    ax[2].set_ylabel('$\sigma$')


def plot_theta2(
        theta2_on, theta2_off, 
        theta2_range=[0, 1.0], n_bins_before_th2_cut=1, 
        theta2_cut=0.04, norm_range=[0.2, 0.5],
        t_elapsed=0 * u.h):
    """
    Plots theta2 distribution for on regions and all off
    regions combined. Li&Ma significance of detection is also
    calculated.

    Parameters
    ----------
    theta2_on: astropy.table.column.Column
    theta2_off: numpy.ndarray
        One column per OFF region
    theta2_range: list
        Range of the horizontal axis in deg^2
    n_bins_before_th2_cut: int
        Binning settings. Instead of the total number of bins,
        N bins up to theta^2 must be provided.
    theta2_cut: float
    norm_range: list
        range of theta^2 from which alpha is calcualted
    t_elapsed: astropy.units.quantity.Quantity

    Returns
    -------

    """

    nbins=round((theta2_range[1]/theta2_cut)*n_bins_before_th2_cut)
    hist_on, bin_edges_on=np.histogram(theta2_on, density=False, bins=nbins, range=theta2_range)
    hist_off, bin_edges_off=np.histogram(theta2_off.flatten(), density=False, bins=nbins, range=theta2_range)
    bin_width=bin_edges_on[1]-bin_edges_off[0]
    bin_center=bin_edges_on[:-1]+(bin_width/2)

    idx_min = (np.abs(bin_edges_on - norm_range[0])).argmin()
    idx_max = (np.abs(bin_edges_on - norm_range[1])).argmin()

    Non_norm = np.sum(hist_on[idx_min:idx_max])
    Noff_norm = np.sum(hist_off[idx_min:idx_max])

    alpha = Non_norm / Noff_norm

    N_on = np.sum(theta2_on<=theta2_cut)
    N_off = np.sum(theta2_off.flatten()<=theta2_cut)
    N_excess = int(N_on - N_off * alpha)

    # number of off regions
    n_off = theta2_off.shape[1]

    stat = WStatCountsStatistic(n_on=N_on, n_off=N_off, alpha=alpha)
    significance_lima = stat.sqrt_ts

    textstr = r'N$_{{\rm on}}$ = {:.0f} '\
                f'\n'\
                r'N$_{{\rm off}}$ = {:.0f} '\
                f'\n'\
                r'N$_{{\rm excess}}$ = {:.0f} '\
                f'\n'\
                r'n$_{{\rm off \, regions}}$ = {:.0f} '\
                f'\n'\
                r'Time = {:.1f}'\
                f'\n'\
                r'LiMa Significance = {:.1f} $\sigma$ '.format(N_on,
                                                          N_off,
                                                          N_excess,
                                                          n_off,
                                                          t_elapsed.to(u.h),
                                                          significance_lima)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.95)
    
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.errorbar(bin_center, hist_on, yerr=np.sqrt(hist_on), fmt='o', label='ON data', ms=10, color='lightcoral')
    ax.errorbar(bin_center, alpha * hist_off, yerr=alpha * np.sqrt(hist_off),fmt='s',label='Background', ms=10, color='cornflowerblue')
    ax.set_xlim(theta2_range[0], theta2_range[1])
    ax.grid(ls='dashed')
    ax.axvline(theta2_cut, color='black',ls='--',alpha=0.75)
    ax.set_xlabel("$\\theta^{2} [deg^{2}]$")
    ax.set_ylabel("Counts")
    ax.legend(bbox_to_anchor=(0.1, 0.95))

    txt = ax.text(0.50, 0.96, textstr, transform=ax.transAxes, fontsize=15,
            verticalalignment='top', bbox=props)


def plot_mc_data(data=None, mc=None, scale=1):
    """
    Plots comparison of MC/data Hillas parameters.

    Parameters
    ----------
    data: pandas.DataFrame with DL1 table
    mc: pandas.DataFrame with DL1 table
    scale: float

    Returns
    -------

    """

    fig, ax = plt.subplots(4, 3, figsize=(20, 25))

    ax[0, 0].hist(data['camera_frame_hillas_intensity'], bins=np.logspace(-1, 5, 100), histtype='step', label='DATA', linewidth=3)
    ax[0, 0].hist(mc['camera_frame_hillas_intensity'], weights=scale*mc["weight"], bins=np.logspace(-1, 5, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[0, 0].set_xscale('log')
    ax[0, 0].set_yscale('log')
    ax[0, 0].set_xlabel('Intensity [p.e.]')
    ax[0, 0].set_ylabel('N')
    ax[0, 0].set_xlim([10**1, 10**5])
    ax[0, 0].grid()
    ax[0, 0].legend()

    ax[0, 1].hist(data['camera_frame_hillas_width'], bins=np.linspace(-0.01, 0.04, 100), histtype='step', label='DATA', linewidth=3)
    ax[0, 1].hist(mc['camera_frame_hillas_width'], weights=scale*mc["weight"], bins=np.linspace(-0.01, 0.04, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[0, 1].set_xlabel('Hillas width')
    ax[0, 1].set_ylabel('N')
    ax[0, 1].set_xlim([-0.01, 0.04])
    ax[0, 1].grid()
    ax[0, 1].legend()

    ax[0, 2].hist(data['camera_frame_hillas_length'], bins=np.linspace(0, 0.08, 100), histtype='step', label='DATA', linewidth=3)
    ax[0, 2].hist(mc['camera_frame_hillas_length'], weights=scale*mc["weight"], bins=np.linspace(0, 0.08, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[0, 2].set_xlabel('Hillas length')
    ax[0, 2].set_ylabel('N')
    ax[0, 2].set_xlim([0, 0.08])
    ax[0, 2].grid()
    ax[0, 2].legend()

    ax[1, 0].hist(data['camera_frame_hillas_skewness'], bins=np.linspace(-2, 2, 100), histtype='step', label='DATA', linewidth=3)
    ax[1, 0].hist(mc['camera_frame_hillas_skewness'], weights=scale*mc["weight"], bins=np.linspace(-2, 2, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[1, 0].set_xlabel('camera_frame_hillas_skewness')
    ax[1, 0].set_ylabel('N')
    ax[1, 0].set_xlim([-2, 2])
    ax[1, 0].grid()
    ax[1, 0].legend()

    ax[1, 1].hist(data['camera_frame_hillas_kurtosis'], bins=np.linspace(0, 5, 100), histtype='step', label='DATA', linewidth=3)
    ax[1, 1].hist(mc['camera_frame_hillas_kurtosis'], weights=scale*mc["weight"], bins=np.linspace(0, 5, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[1, 1].set_xlabel('camera_frame_hillas_kurtosis')
    ax[1, 1].set_ylabel('N')
    ax[1, 1].set_xlim([0, 5])
    ax[1, 1].grid()
    ax[1, 1].legend()

    ax[1, 2].hist(data['camera_frame_timing_slope'], bins=np.linspace(-150, 150, 100), histtype='step', label='DATA', linewidth=3)
    ax[1, 2].hist(mc['camera_frame_timing_slope'], weights=scale*mc["weight"], bins=np.linspace(-150, 150, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[1, 2].set_xlabel('Time gradient')
    ax[1, 2].set_ylabel('N')
    ax[1, 2].set_xlim([-150, 150])
    ax[1, 2].grid()
    ax[1, 2].legend()

    ax[2, 0].hist(data['camera_frame_hillas_x'], bins=np.linspace(-0.75, 0.75, 100), histtype='step', label='DATA', linewidth=3)
    ax[2, 0].hist(mc['camera_frame_hillas_x'], weights=scale*mc["weight"], bins=np.linspace(-0.75, 0.75, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[2, 0].set_xlabel('cog x')
    ax[2, 0].set_ylabel('N')
    ax[2, 0].set_xlim([-0.75, 0.75])
    ax[2, 0].grid()
    ax[2, 0].legend()

    ax[2, 1].hist(data['camera_frame_hillas_y'], bins=np.linspace(-0.75, 0.75, 100), histtype='step', label='DATA', linewidth=3)
    ax[2, 1].hist(mc['camera_frame_hillas_y'], weights=scale*mc["weight"], bins=np.linspace(-0.75, 0.75, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[2, 1].set_xlabel('cog y')
    ax[2, 1].set_ylabel('N')
    ax[2, 1].set_xlim([-0.75, 0.75])
    ax[2, 1].grid()
    ax[2, 1].legend()

    ax[2, 2].hist(data['camera_frame_hillas_phi'], bins=np.linspace(-190, 190, 100), histtype='step', label='DATA', linewidth=3)
    ax[2, 2].hist(mc['camera_frame_hillas_phi'], weights=scale*mc["weight"], bins=np.linspace(-190, 190, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[2, 2].set_xlabel('phi')
    ax[2, 2].set_ylabel('N')
    ax[2, 2].set_xlim([-190, 190])
    ax[2, 2].grid()
    ax[2, 2].legend()

    ax[3, 0].hist(data['camera_frame_hillas_psi'], bins=np.linspace(-250, 250, 100), histtype='step', label='DATA', linewidth=3)
    ax[3, 0].hist(mc['camera_frame_hillas_psi'], weights=scale*mc["weight"], bins=np.linspace(-250, 250, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[3, 0].set_xlabel('psi')
    ax[3, 0].set_ylabel('N')
    ax[3, 0].set_xlim([-250, 250])
    ax[3, 0].grid()
    ax[3, 0].legend()

    ax[3, 1].hist(data['camera_frame_hillas_r'], bins=np.linspace(0, 0.6, 100), histtype='step', label='DATA', linewidth=3)
    ax[3, 1].hist(mc['camera_frame_hillas_r'], weights=scale*mc["weight"], bins=np.linspace(0, 0.6, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[3, 1].set_xlabel('r')
    ax[3, 1].set_ylabel('N')
    ax[3, 1].set_xlim([0, 0.6])
    ax[3, 1].grid()
    ax[3, 1].legend()

    ax[3, 2].hist(data['leakage_intensity_width_2'], bins=np.linspace(0, 1, 100), histtype='step', label='DATA', linewidth=3)
    ax[3, 2].hist(mc['leakage_intensity_width_2'], weights=scale*mc["weight"], bins=np.linspace(0, 1, 100), histtype='step', label='MC (scaled)', linewidth=3)
    ax[3, 2].set_xlabel('leakage_intensity_width_2')
    ax[3, 2].set_ylabel('N')
    ax[3, 2].set_xlim([0, 1])
    plt.yscale('log')
    ax[3, 2].grid()
    ax[3, 2].legend()

    fig.tight_layout()


def plot_preliminary(ax=None, fontsize=60):
    ax.text(0.5, 0.5, 'PRELIMINARY', transform=ax.transAxes,
        fontsize=fontsize, color='grey', alpha=0.3,
        ha='center', va='center', rotation=30)


def plot_energy_resolution(
        e_bins, e_tables=None, labels=None, 
        markers=None, axes=None, plot=None, 
        preliminary=False,
        skip_bins_first=0, skip_bins_last=0,
        titles=['Energy resolution', 'Energy bias']):
    """
    Plots energy resolution and bias.

    Parameters
    ----------
    e_bins: astropy.table.table.QTable
    e_tables: list of astropy.table.table.Table
    labels: list of strings
    markers: list of strings
        List of markers 
    axes: numpy.ndarray of matplotlib.axes._axes.Axes
        Axes to plot the figure in 
    plot: str
        If not provided with two axes, user must
        select whether to plot \'resolution\'
        or \'bias\'
    preliminary: bool
        Add preliminary watermark
    skip_bins_first: float
        Skip certain number of first bins, e.g. with very low number of events
    skip_bins_last: float
        Skip certain number of last bins
    titles: list of strings
        Plot titles

    Returns
    -------

    """

    energy_center = ctaplot.ana.logbin_mean(e_bins['energy_bins'])
    
    marker_size = 6
    if markers == None:
        markers = len(labels) * ['o']
        
    if len(axes) == 2:
        
        if preliminary:
            plot_preliminary(ax=axes[0])
            plot_preliminary(ax=axes[1])

        for table, label, marker in zip(e_tables, labels, markers):

            energy_center_plot = energy_center[skip_bins_first:len(energy_center)-skip_bins_last]
            table_plot = table[skip_bins_first:len(energy_center)-skip_bins_last]

            axes[0].errorbar(energy_center_plot, table_plot['energy_res'],
                    xerr=(energy_center_plot - e_bins['energy_bins'][skip_bins_first:len(e_bins)-1-skip_bins_last], e_bins['energy_bins'][skip_bins_first+1:len(e_bins)-skip_bins_last] - energy_center_plot),
                    yerr=(table_plot['energy_res'] - table_plot['energy_res_err_lo'], table_plot['energy_res_err_hi'] - table_plot['energy_res']),
                    fmt=marker,
                    label=label,
                    markersize=marker_size
                    )
            axes[1].errorbar(energy_center_plot, 
                    table_plot['energy_bias'], 
                    xerr=(energy_center_plot - e_bins['energy_bins'][skip_bins_first:len(e_bins)-1-skip_bins_last], e_bins['energy_bins'][skip_bins_first+1:len(e_bins)-skip_bins_last] - energy_center_plot),
                    fmt=marker, 
                    label=label
                    )

        energy_center_unit = energy_center.unit.to_string("latex")
        axes[0].legend()
        axes[1].legend()
        axes[0].set_ylabel(r"$(\Delta E/E_\mathrm{true})_{68}$")
        axes[0].set_xlabel(r'$E_\mathrm{true}$' + ' [' + energy_center_unit + ']')
        axes[0].set_xscale('log')
        if len(titles):
            axes[0].set_title(titles[0])
            axes[1].set_title(titles[1])
        axes[0].set_xlim([min(e_bins['energy_bins'].value), max(e_bins['energy_bins'].value)])
        axes[0].set_ylim([0, 0.7])
        axes[0].grid(True, which='both')
        axes[0].label_outer()
        axes[1].set_xscale('log')            
        axes[1].grid(True, which='both')
        axes[1].set_ylabel("Energy bias")
        axes[1].set_xlim([min(e_bins['energy_bins'].value), max(e_bins['energy_bins'].value)])
        axes[1].set_ylim([-0.3, 0.3])
        axes[1].set_ylabel(r"bias (median($E_\mathrm{reco}/E_\mathrm{true}$ - 1)")
        axes[1].set_xlabel(r'$E_\mathrm{true}$' + ' [' + energy_center_unit + ']')
    else:
        if plot == 'resolution':
            
            if preliminary:
                plot_preliminary(ax=axes[0])
            
            for table, label, marker in zip(e_tables, labels, markers):

                energy_center_plot = energy_center[skip_bins_first:len(energy_center)-skip_bins_last]
                table_plot = table[skip_bins_first:len(energy_center)-skip_bins_last]

                axes[0].errorbar(energy_center_plot, table_plot['energy_res'],
                    xerr=(energy_center_plot - e_bins['energy_bins'][skip_bins_first:len(e_bins)-1-skip_bins_last], e_bins['energy_bins'][skip_bins_first+1:len(e_bins)-skip_bins_last] - energy_center_plot),
                    yerr=(table_plot['energy_res'] - table_plot['energy_res_err_lo'], table_plot['energy_res_err_hi'] - table_plot['energy_res']),
                    fmt=marker,
                    label=label,
                    markersize=marker_size
                    )

            energy_center_unit = energy_center.unit.to_string("latex")
            axes[0].legend()
            axes[0].set_ylabel(r"$(\Delta E/E_\mathrm{true})_{68}$")
            axes[0].set_xlabel(r'$E_\mathrm{true}$' + ' [' + energy_center_unit + ']')
            axes[0].set_xscale('log')
            if len(titles):
                axes[0].set_title(titles[0])
            axes[0].set_xlim([min(e_bins['energy_bins'].value), max(e_bins['energy_bins'].value)])
            axes[0].set_ylim([0, 0.8])
            axes[0].grid(True, which='both')
            axes[0].label_outer()
        elif plot == 'bias':
            
            if preliminary:
                plot_preliminary(ax=axes[0])

            for table, label, marker in zip(e_tables, labels, markers):

                energy_center_plot = energy_center[skip_bins_first:len(energy_center)-skip_bins_last]
                table_plot = table[skip_bins_first:len(energy_center)-skip_bins_last]

                axes[0].errorbar(energy_center_plot, 
                                 table_plot['energy_bias'], 
                                 xerr=(energy_center_plot - e_bins['energy_bins'][skip_bins_first:len(e_bins)-1-skip_bins_last], e_bins['energy_bins'][skip_bins_first+1:len(e_bins)-skip_bins_last] - energy_center_plot),
                                 fmt=marker, 
                                 label=label,
                                  markersize=marker_size
                                )

            energy_center_unit = energy_center.unit.to_string("latex")
            axes[0].legend()
            axes[0].set_xscale('log')
            if len(titles):
                axes[0].set_title(titles[0])
            axes[0].grid(True, which='both')
            axes[0].set_ylabel("Energy bias")
            axes[0].set_xlim([min(e_bins['energy_bins'].value), max(e_bins['energy_bins'].value)])
            axes[0].set_ylim([-0.3, 0.3])
            axes[0].set_ylabel(r"bias (median($E_\mathrm{reco}/E_\mathrm{true}$ - 1)")
            axes[0].set_xlabel(r'$E_\mathrm{true}$' + ' [' + energy_center_unit + ']')

        else: print('Set plot=resolution or plot=bias.')


def plot_angular_resolution(
        e_bins, a_tables=None, labels=None, 
        markers=None, ax=None, preliminary=False,
        skip_bins_first=0, skip_bins_last=0,
        title='Angular resolution'):
    """
    Plots energy resolution and bias.

    Parameters
    ----------
    e_bins: astropy.table.table.QTable
    a_tables: list of astropy.table.table.Table
    labels: list of strings
    markers: list of strings
        List of markers 
    ax: matplotlib.axes._axes.Axes
        Axis to plot the figure in 
    preliminary: bool
        Add preliminary watermark
    skip_bins_first: float
        Skip certain number of first bins, e.g. with very low number of events
    skip_bins_last: float
        Skip certain number of last bins
    title: string
        Plot title

    Returns
    -------

    """

    energy_center = ctaplot.ana.logbin_mean(e_bins['energy_bins'])
    
    marker_size = 6
    if markers == None:
        markers = len(labels) * ['o']

    if preliminary:
        plot_preliminary(ax=ax)

    for table, label, marker in zip(a_tables, labels, markers):
        
        energy_center_plot = energy_center[skip_bins_first:len(energy_center)-skip_bins_last]
        table_plot = table[skip_bins_first:len(energy_center)-skip_bins_last]
        
        ax.errorbar(energy_center_plot, 
                    table_plot['angular_res'], 
                    xerr=(energy_center_plot - e_bins['energy_bins'][skip_bins_first:len(e_bins)-1-skip_bins_last], e_bins['energy_bins'][skip_bins_first+1:len(e_bins)-skip_bins_last] - energy_center_plot),
                    yerr=(table_plot['angular_res'] - table_plot['angular_res_err_lo'], table_plot['angular_res_err_hi'] - table_plot['angular_res']),
                    fmt=marker, 
                    markersize=marker_size,
                    label=label)

    energy_center_unit = energy_center.unit.to_string("latex")
    ax.set_ylabel('Angular Resolution [deg]')
    ax.set_xlabel(r'$E_\mathrm{true}$' + ' ['+ energy_center_unit + ']')
    ax.set_xscale('log')
    if len(title):
        ax.set_title(title)
    ax.grid(True, which='both')
    ax.set_xlim([min(e_bins['energy_bins'].value), max(e_bins['energy_bins'].value)])
    ax.set_ylim([0, 0.6])
    ax.legend()


def plot_roc(
        roc_tables=None, labels=None, ax=None, 
        gammaness_cut=None, preliminary=False, 
        linestyles=None, title='ROC curve'):
    """
    Plots energy resolution and bias.

    Parameters
    ----------
    roc_tables: list of astropy.table.table.Table
    labels: list of strings
    ax: matplotlib.axes._axes.Axes
        Axis to plot the figure in 
    gammaness_cut: float
        Highlights one point on the ROC curve
    preliminary: bool
        Add preliminary watermark
    linestyles:
        List of linestyles
    title: string
        Plot title

    Returns
    -------

    """

    if preliminary:
        plot_preliminary(ax=ax)

    if linestyles == None:
        linestyles = len(labels) * ['-']

    for table, label, linestyle in zip(roc_tables, labels, linestyles):
        
        ax.plot(table['false_positive_rate'], 
                table['true_positive_rate'], 
                label=label + ', AUC=' + str(get_auc(table)),
                linestyle=linestyle
               )
        if gammaness_cut is not None:
            ax.plot(table['false_positive_rate'][np.argmin(abs(table['threshold'] - gammaness_cut))], 
                    table['true_positive_rate'][np.argmin(abs(table['threshold'] - gammaness_cut))],
                    'k.')
    ax.set_ylabel('True Positive Rate')
    ax.set_xlabel('False Positive Rate')
    if len(title):
        ax.set_title(title)
    ax.grid(True, which='both')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.legend()


def plot_sensitivity(
        sens_tables=None, labels=None, markers=None, 
        ax=None, xlim=[10**-1, 10**3], ylim=[10**-12, 10**-10], 
        bands=False, yerrors=False, preliminary=False,
        legend_fontsize=10):
    """
    Plots energy resolution and bias.

    Parameters
    ----------
    sens_tables: list of astropy.table.table.Table
    labels: list of strings
    markers:
        List of markers
    ax: matplotlib.axes._axes.Axes
        Axis to plot the figure in 
    xlim: list of floats
    ylim: list of floats
    bands: bool
        Plot uncertainties as bands
    yerrors: bool
        Plot uncertainties
    preliminary: bool
        Add preliminary watermark
    legend_fontsize: int
        Legend font size

    Returns
    -------

    """

    if preliminary:
        plot_preliminary(ax=ax)

    marker_size = 6
    if markers == None:
        markers = len(labels) * ['o']
    
    for table, label, marker in zip(sens_tables, labels, markers):

        if bands:
            ax.fill_between(table['energy'],
                            table['flux_sensitivity']-table['flux_sensitivity_err_minus'],
                            table['flux_sensitivity']+table['flux_sensitivity_err_plus'],
                            alpha=0.5, 
                            antialiased=True,
                            label=label
                           )
        else:
            if yerrors:
                ax.errorbar(table['energy'], table['flux_sensitivity'],
                            xerr=(table['energy'] - table['energy_min'], table['energy_max'] - table['energy']),
                            yerr=(table['flux_sensitivity_err_minus'], table['flux_sensitivity_err_plus']),
                            fmt=marker,
                            markersize=marker_size,
                            label=label
                            )
            else:
                ax.errorbar(table['energy'], table['flux_sensitivity'],
                            xerr=(table['energy'] - table['energy_min'], table['energy_max'] - table['energy']),
                            fmt=marker,
                            markersize=marker_size,
                            label=label
                            )
    units_flux = table['flux_sensitivity'].unit
    units_energy = table['energy'].unit.to_string("latex")
    ax.set_ylabel(rf'Flux sensitivity [{units_flux.to_string("latex")}]')
    ax.set_xlabel(r'$E_\mathrm{reco}$ '+ '[' + units_energy + ']')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, which='both')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.legend(fontsize=legend_fontsize)


def get_auc(table):
    """
    Calculates the \'Area under the curve\'

    Parameters
    ----------
    table: list of astropy.table.table.Table

    Returns
    -------
    auc: float

    """

    d_fpr = table['false_positive_rate'][1:]-table['false_positive_rate'][:-1]
    auc = round(sum(table['true_positive_rate'][:-1] * d_fpr), 3)
    return auc


def plot_source_sed(
        spectrum=None, ax=None, fraction=None, 
        linestyle="-", color='grey', 
        flux_units=u.TeV / (u.cm ** 2 * u.s), 
        alpha=None):
    """
    Plot SED of given gamma-ray source.

    Parameters
    ----------
    spectrum: string
        Designation of source spectra as defined in
        sst1mpipe.performance.spectra
    ax: matplotlib.axes._axes.Axes
        Axis to plot the figure in 
    fraction: float
        Scaling of plotted SED
    linestyle: string
    color: string
    flux_units: astropy.units.CompositeUnit
    alpha: float

    Returns
    -------
    h: list of matplotlib.lines.Line2D

    """
    
    target_spectrum = globals()[spectrum]
    energy_smooth = np.logspace(-1, 3, 200) * u.TeV
    if fraction is not None:
        flux = fraction*target_spectrum(energy_smooth)
        h = ax.plot(energy_smooth, (flux * energy_smooth**2).to(flux_units), color=color, label=str(fraction)+" "+spectrum, linestyle=linestyle, alpha=alpha)
    else:
        flux = target_spectrum(energy_smooth)
        h = ax.plot(energy_smooth, (flux * energy_smooth**2).to(flux_units), color=color, label=spectrum, linestyle=linestyle, alpha=alpha)
    return h


def plot_lhaaso_sens(ax=None):
    """
    Plot LHAASO sensitivity from Fig 2 in
    LHAASO col 2018, https://link.springer.com/article/10.1007/s41605-018-0037-3

    Parameters
    ----------
    ax: matplotlib.axes._axes.Axes
        Axis to plot the figure in 

    Returns
    -------

    """
    
    sens = np.array(
        [
            [153.6634615886556, 3.754698473106022e-11],
            [240.6985340697518, 1.2035478358250054e-11],
            [489.66566456144955, 3.815576123757241e-12],
            [895.8580662767178, 1.4363464415763194e-12],
            [1548.3791387970177, 1.0744385479437988e-12],
            [2697.7334209590035, 5.76317661191093e-13],
            [4587.421608158554, 4.956365252550058e-13],
            [8746.048036076914, 4.4013376413891905e-13],
            [15878.079331542003, 4.542143822336942e-13],
            [23684.84614054241, 1.9452468317768614e-13],
            [41605.06043884821, 1.0775285670493913e-13],
            [74285.6092286208, 5.596530590109233e-14],
            [132670.00037320558, 5.1327626820170326e-14],
            [233110.56709684362, 5.1293938735424906e-14],
            [419813.2944919163, 7.07207938612841e-14],
            [750173.7128410875, 2.2273960279095398e-13]
        ]
    )
    
    energy = sens[:, 0] * u.GeV
    diff_sens = sens[:, 1] * u.erg / (u.cm ** 2 * u.s)
    ax.plot(energy.to(u.TeV), diff_sens.to(u.TeV / (u.cm ** 2 * u.s)), label='LHAASO (1yr)')


def plot_hawc_sens(ax=None):
    """
    Plot HAWC sensitivity from Fig 2 in
    HAWC col 2016, https://www.sciencedirect.com/science/article/pii/S2405601415005295?via%3Dihub

    Parameters
    ----------
    ax: matplotlib.axes._axes.Axes
        Axis to plot the figure in 

    Returns
    -------

    """

    sens = np.array(
        [
            [312.70754082889096, 2.4821784185383784e-11],
            [406.0795379904598, 1.8875085975553973e-11],
            [506.2332582968051, 1.474394701077434e-11],
            [620.8646990993665, 1.1641424773526121e-11],
            [725.0424589638405, 9.492949409859737e-12],
            [846.6966588026186, 7.65837996385966e-12],
            [1005.0450429188891, 6.1123061818792504e-12],
            [1212.6528159820657, 4.826204726479716e-12],
            [1511.765088843364, 3.935208174578947e-12],
            [1931.3762630218828, 3.073832069878877e-12],
            [2487.7423643653283, 2.506256860689843e-12],
            [3310.844065311097, 2.0877239084498234e-12],
            [4299.773424925843, 1.894991404163835e-12],
            [5769.647682799492, 1.7762445377550779e-12],
            [7554.718249174407, 1.7379921555771015e-12],
            [10055.346794570081, 1.833156795144318e-12],
            [13274.931754220937, 1.975486747521003e-12],
            [16962.4764589694, 2.2704967264869583e-12],
            [21498.331107678358, 2.6949433051918457e-12],
            [27025.298571268788, 3.1646321190090285e-12],
            [34533.27702513294, 3.837659388955406e-12],
            [42361.65637211649, 4.806269778115397e-12],
            [51121.8468366289, 5.828814924606857e-12],
            [62709.80880656513, 7.068774421636179e-12],
            [76924.45331881678, 8.572509587332269e-12],
            [103236.73479173059, 1.1326629017510099e-11]
        ]
    )
    
    energy = sens[:, 0] * u.GeV
    diff_sens = sens[:, 1] * u.erg / (u.cm ** 2 * u.s)
    ax.plot(energy.to(u.TeV), diff_sens.to(u.TeV / (u.cm ** 2 * u.s)), label='HAWC (1yr)')


def plot_veritas_sens(ax=None):
    """
    Plot VERITAS sensitivity from Fig 22 in
    https://ui.adsabs.harvard.edu/abs/2020NCimR..43..281C/abstract

    Parameters
    ----------
    ax: matplotlib.axes._axes.Axes
        Axis to plot the figure in 

    Returns
    -------

    """

    sens = np.array(
        [[0.13445071561209565, 2.7989765487454044e-12],
        [0.17295336993434499, 1.8696066912593394e-12],
        [0.1935506561970808, 1.6434612156435331e-12],
        [0.2612780744931726, 1.1740721737263156e-12],
        [0.2939649289555213, 1.0613730003060555e-12],
        [0.43234932738121645, 8.29297135071508e-13],
        [0.48643782894575727, 7.66682880180691e-13],
        [0.7348543677167647, 6.16050434416116e-13],
        [0.8446982037963724, 6.057315536295814e-13],
        [1.3319593877373204, 6.368801273889653e-13],
        [1.5642242378556541, 6.512497889588396e-13],
        [2.4665414804650796, 6.9245455429238e-13],
        [2.8200822827770358, 7.080896571438344e-13],
        [4.375930237948393, 8.143281310482352e-13],
        [5.00315256501573, 8.956257737944712e-13],
        [7.358393581016099, 1.1521624697566649e-12],
        [8.323433257799676, 1.267195729052601e-12],
        [12.307467645029151, 1.5762710035213683e-12]]
        )

    energy = sens[:, 0] * u.TeV
    sens = sens[:, 1] * u.erg * u.cm**-2 * u.s**-1
    ax.plot(energy.to(u.TeV), sens.to(u.TeV / (u.cm ** 2 * u.s)), label='VERITAS (50h)')


def plot_hess_sens(ax=None):
    """
    Plot H.E.S.S. sensitivity from Fig 22 in
    https://ui.adsabs.harvard.edu/abs/2020NCimR..43..281C/abstract

    Parameters
    ----------
    ax: matplotlib.axes._axes.Axes
        Axis to plot the figure in 

    Returns
    -------

    """

    sens = np.array(
        [
            [0.20568922865735784, 2.6463988455392184e-12],
            [0.2275964287543443, 1.999183315947727e-12],
            [0.2518631968873063, 1.5795681488794446e-12],
            [0.2847625062147546, 1.255047521897012e-12],
            [0.345198305137716, 9.9719931890071e-13],
            [0.4141274799274541, 9.065026168181143e-13],
            [0.5049282682962686, 8.618755879418285e-13],
            [0.6254886560364484, 7.491059385323583e-13],
            [0.800195192466779, 6.621403604501413e-13],
            [1.0573844506855754, 6.402266371250526e-13],
            [1.3749261230006593, 6.155755533998118e-13],
            [1.787781601804662, 5.852708333830761e-13],
            [2.168650543699582, 6.190381546984956e-13],
            [2.520548717686202, 6.925298375355217e-13],
            [2.9613214549236844, 7.967825779923441e-13],
            [3.5732747234716133, 8.715989270799347e-13],
            [4.747632044258728, 8.814319610804054e-13],
            [5.979922560466048, 9.750744098492465e-13],
            [7.413055700272808, 1.1537703301130702e-12]
        ]
    )

    energy = sens[:, 0] * u.TeV
    sens = sens[:, 1] * u.erg * u.cm**-2 * u.s**-1
    ax.plot(energy.to(u.TeV), sens.to(u.TeV / (u.cm ** 2 * u.s)), label='H.E.S.S. (50h)')




def plot_astri_sens(ax=None):
    """
    Plot ASTRI sensitivity from Fig 4 in
    ASTRI col 2022, https://arxiv.org/pdf/2208.03177.pdf

    Parameters
    ----------
    ax: matplotlib.axes._axes.Axes
        Axis to plot the figure in 

    Returns
    -------

    """

    sens = np.array(
        [
            [0.4004203818792009, 6.844587889113004e-11],
            [0.6314244219820897, 3.83378946293183e-12],
            [0.9952911474633526, 1.8684461844733356e-12],
            [1.5798789597490586, 1.005110850202892e-12],
            [2.4937016271797585, 7.812581684540606e-13],
            [3.9381551209028465, 7.266873409566725e-13],
            [6.261739469937325, 6.943727773752925e-13],
            [9.893200990187538, 7.523540469549305e-13],
            [15.738620533172897, 8.602798812796036e-13],
            [25.042413157337197, 1.0474845267985117e-12],
            [39.85131003247597, 1.3339777003835077e-12],
            [63.002477617123105, 1.7928544826328416e-12],
            [99.58739170774423, 2.283230666974204e-12],
            [158.55344222243886, 3.417636857094684e-12]
        ]
    )
    energy = sens[:, 0] * u.TeV
    diff_sens = sens[:, 1] * u.erg / (u.cm ** 2 * u.s)
    ax.plot(energy.to(u.TeV), diff_sens.to(u.TeV / (u.cm ** 2 * u.s)), label='ASTRI Mini-Array (50h)')


def plot_theta2_dl3(ax=None, theta2_axis=None, counts_on=None, counts_off=None, alpha=None, theta_cut=None, event_counts=None):

    ax.errorbar(theta2_axis.center, counts_on, yerr=np.sqrt(counts_on), fmt='o', ms=8, label='ON data')
    ax.errorbar(theta2_axis.center, alpha*counts_off, yerr=alpha*np.sqrt(counts_off), fmt='o', ms=8, label='OFF data')
    ax.set_xlabel("$\\theta^{2} [deg^{2}]$")
    ax.set_ylabel("Counts")
    ax.grid(ls='dashed')
    ax.axvline(theta_cut.to_value()**2, color='black',ls='--',alpha=0.75)
    ax.set_xlim(theta2_axis.bounds[0].value, theta2_axis.bounds[1].value)

    textstr = r'N$_{{\rm on}}$ = {:.0f} '\
                f'\n'\
                r'N$_{{\rm off}}$ = {:.0f} '\
                f'\n'\
                r'N$_{{\rm excess}}$ = {:.0f} '\
                f'\n'\
                r'n$_{{\rm off \, regions}}$ = {:.0f} '\
                f'\n'\
                r'Time = {:.1f}'\
                f'\n'\
                r'LiMa Significance = {:.1f} $\sigma$ '.format(event_counts.N_on,
                                                        event_counts.N_off,
                                                        event_counts.N_excess,
                                                        event_counts.n_off_regions,
                                                        event_counts.t_elapsed,
                                                        event_counts.significance_lima)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.95)
    txt = ax.text(0.50, 0.96, textstr, transform=ax.transAxes, fontsize=13,
                verticalalignment='top', bbox=props)
    ax.legend(bbox_to_anchor=(0.1, 0.95), fontsize=10)