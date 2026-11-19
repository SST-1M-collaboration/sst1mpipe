import itertools

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve
from astropy.convolution.kernels import Gaussian2DKernel
from astropy.visualization.wcsaxes import SphericalCircle
from gammapy.stats import WStatCountsStatistic


def plot_count_maps(
        data, bins_raw=20, bins_conv=100,
        conv_kernel_deg=0.2, range_deg=4,
        source=None, theta2_cut=None,
        wobbles=None, off_radec_w=None):
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

    if wobbles is None:
        wobbles = []
    if off_radec_w is None:
        off_radec_w = []

    range_d = [
                       [source.ra.value-range_deg/2/np.cos(np.deg2rad(source.dec.value)),
                        source.ra.value+range_deg/2/np.cos(np.deg2rad(source.dec.value))],
                       [source.dec.value-range_deg/2,
                        source.dec.value+range_deg/2]
                   ]

    fig, ax = plt.subplots(1, 2, figsize=(17,7))

    ax[0].hist2d(data.reco_ra,
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
        ax[1].pcolor(xedges[edge:-edge], yedges[edge:-edge], image[edge:-edge, edge:-edge].T)
    else:
        ax[1].pcolor(xedges, yedges, image.T)

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
    yield from itertools.cycle(iterable)


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
    ax[0].set_ylabel(r'$\sigma$')

    ax[1].plot(times, alphas, 'b.')
    ax[1].grid()
    ax[1].set_xlabel(r'$t \, [s]$')
    ax[1].set_ylabel(r'$\alpha$')

    times = np.array(times)
    sigma = np.array(sigma)
    obsid = np.array(obsid)

    colors = ['r', 'b', 'g', 'y', 'm']

    for i, c in enumerate(cycle(colors)):
        if i < len(nights)-1:
            mask = (obsid >= nights[i]) & (obsid <= nights[i+1])
        else:
            mask = (obsid >= nights[i])
        ax[2].plot(times[mask]/3600, sigma[mask], '.', color=c)

        ax[2].text(times[mask][0]/3600, sigma[mask][0]-1, str(nights[i])[4:8], fontsize=12, color=c, fontweight='bold')

    ax[2].plot(times/3600, (max(sigma) / np.sqrt(times[-1]/3600)) * np.sqrt(times/3600), 'y-')
    ax[2].grid()
    ax[2].set_xlabel('t [h]')
    ax[2].set_ylabel(r'$\sigma$')


def plot_theta2(
        theta2_on, theta2_off,
        theta2_range=None, n_bins_before_th2_cut=1,
        theta2_cut=0.04, norm_range=None,
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

    if theta2_range is None:
        theta2_range = [0, 1.0]
    if norm_range is None:
        norm_range = [0.2, 0.5]

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

    textstr = rf'N$_{{\rm on}}$ = {N_on:.0f} '\
                '\n'\
                rf'N$_{{\rm off}}$ = {N_off:.0f} '\
                '\n'\
                rf'N$_{{\rm excess}}$ = {N_excess:.0f} '\
                '\n'\
                rf'n$_{{\rm off \, regions}}$ = {n_off:.0f} '\
                '\n'\
                rf'Time = {t_elapsed.to(u.h):.1f}'\
                '\n'\
                rf'LiMa Significance = {significance_lima:.1f} $\sigma$ '

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

    ax.text(0.50, 0.96, textstr, transform=ax.transAxes, fontsize=15,
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


def plot_theta2_dl3(ax=None, theta2_axis=None, counts_on=None, counts_off=None, alpha=None, theta_cut=None, event_counts=None):

    ax.errorbar(theta2_axis.center, counts_on, yerr=np.sqrt(counts_on), fmt='o', ms=8, label='ON data')
    ax.errorbar(theta2_axis.center, alpha*counts_off, yerr=alpha*np.sqrt(counts_off), fmt='o', ms=8, label='OFF data')
    ax.set_xlabel("$\\theta^{2} [deg^{2}]$")
    ax.set_ylabel("Counts")
    ax.grid(ls='dashed')
    ax.axvline(theta_cut.to_value()**2, color='black',ls='--',alpha=0.75)
    ax.set_xlim(theta2_axis.bounds[0].value, theta2_axis.bounds[1].value)

    textstr = rf'N$_{{\rm on}}$ = {event_counts.N_on:.0f} '\
                '\n'\
                rf'N$_{{\rm off}}$ = {event_counts.N_off:.0f} '\
                '\n'\
                rf'N$_{{\rm excess}}$ = {event_counts.N_excess:.0f} '\
                '\n'\
                rf'n$_{{\rm off \, regions}}$ = {event_counts.n_off_regions:.0f} '\
                '\n'\
                rf'Time = {event_counts.t_elapsed:.1f}'\
                '\n'\
                rf'LiMa Significance = {event_counts.significance_lima:.1f} $\sigma$ '

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.95)
    ax.text(0.50, 0.96, textstr, transform=ax.transAxes, fontsize=13,
                verticalalignment='top', bbox=props)
    ax.legend(bbox_to_anchor=(0.1, 0.95), fontsize=10)
