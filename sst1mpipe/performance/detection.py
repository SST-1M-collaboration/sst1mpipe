"""
Time to detect a point source with the point-like IRFs of ``ctapipe-compute-irf``
(``irf_point_like.fits.gz``): the rates of the signal (spectrum of the source folded with the effective
area and energy dispersion, the effective area including the direction cut) and of the background
(``BACKGROUND`` in the ON region of radius ``RAD_MAX``) above a minimum reconstructed energy, and the
time to reach 5 sigma (WStat, Li & Ma) with ``n_off`` OFF regions.

The counts of the source and of the background grow linearly with the time: so does the test statistic
of Li & Ma for the expected counts, and the time to reach ``sigma`` is ``sigma² / TS(1 s)``.
"""
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.table import QTable
from gammapy.irf import Background2D, EffectiveAreaTable2D, EnergyDispersion2D, RadMax2D
from gammapy.maps import MapAxis
from gammapy.stats import WStatCountsStatistic

from sst1mpipe.performance import spectra

__all__ = ["DetectionRates", "detection_rates", "source_spectrum", "time_to_detection", "time_to_detection_table"]


def source_spectrum(name):
    """Spectrum of a source of `sst1mpipe.performance.spectra` (e.g. ``CRAB_HEGRA``)"""
    if name not in spectra.__all__ or "UNIT" in name:
        raise ValueError(f"Unknown source {name}: one of {[n for n in spectra.__all__ if 'UNIT' not in n]}")
    return getattr(spectra, name)


class DetectionRates:
    """
    Rates of the signal and of the background in the ON region per reconstructed energy bin

    Attributes
    ----------
    energy_edges: Quantity (TeV), reconstructed energy bins
    signal, background: Quantity (s-1), per bin
    """

    def __init__(self, energy_edges, signal, background):
        self.energy_edges = energy_edges
        self.signal = signal
        self.background = background

    def above(self, energy_min):
        """Signal and background rates (s-1) of the bins above ``energy_min`` (bins with their low edge above)"""
        mask = self.energy_edges[:-1] >= u.Quantity(energy_min, u.TeV) * (1 - 1e-6)
        return self.signal[mask].sum(), self.background[mask].sum()


def detection_rates(irf_file, spectrum, offset=None):
    """
    Rates of the signal of a point source of ``spectrum`` (e.g. `source_spectrum`) and of the background
    per bin of reconstructed energy (bins of the ``BACKGROUND`` HDU) with the point-like IRFs ``irf_file``,
    at the offset ``offset`` from the center of the field of view (default: center of the first offset bin)
    """
    irf_file = Path(irf_file)
    aeff = EffectiveAreaTable2D.read(irf_file, hdu="EFFECTIVE AREA")
    edisp = EnergyDispersion2D.read(irf_file, hdu="ENERGY DISPERSION")
    rad_max = RadMax2D.read(irf_file, hdu="RAD_MAX")
    background = Background2D.read(irf_file, hdu="BACKGROUND")
    if offset is None:
        offset = aeff.axes["offset"].center[0]

    energy_true = aeff.axes["energy_true"]
    energy = background.axes["energy"]
    # signal: flux per true energy bin x effective area (direction cut included), migrated to reconstructed energy
    flux = _integrate(spectrum, energy_true.edges)
    rate_true = (flux * aeff.evaluate(energy_true=energy_true.center, offset=offset)).to(u.s**-1)
    kernel = edisp.to_edisp_kernel(offset=offset, energy_axis_true=energy_true,
                                   energy_axis=MapAxis.from_edges(energy.edges, name="energy"))
    signal = (rate_true[:, np.newaxis] * kernel.pdf_matrix).sum(axis=0)

    # background: rate per solid angle and energy x energy bin width x solid angle of the ON region
    radius = rad_max.evaluate(energy=energy.center, offset=offset)
    solid_angle = 2 * np.pi * (1 - np.cos(radius.to_value(u.rad))) * u.sr
    rate_background = (background.evaluate(energy=energy.center, offset=offset) * energy.bin_width * solid_angle).to(u.s**-1)
    return DetectionRates(energy.edges.to(u.TeV), signal.to(u.s**-1), rate_background)


def _integrate(spectrum, edges, n=20):
    """Integral of dN/dE in each bin of ``edges`` (log trapezoids)"""
    integrals = []
    for low, high in zip(edges[:-1], edges[1:], strict=True):
        energy = np.geomspace(low, high, n)
        integrals.append(np.trapezoid(spectrum(energy), energy))
    return u.Quantity(integrals)


def time_to_detection(signal_rate, background_rate, n_off=5, sigma=5, min_signal_events=10,
                      min_excess_over_background=0.05):
    """
    Time to detect a source of signal rate ``signal_rate`` with the background rate ``background_rate``
    in the ON region (``n_off`` OFF regions of the same size), with the conditions of the sensitivity of
    ``ctapipe-compute-irf`` (pyirf): test statistic of Li & Ma (WStat) of the expected counts reaching
    ``sigma²``, at least ``min_signal_events`` signal events, and a signal at least
    ``min_excess_over_background`` times the background (else: infinite)
    """
    signal_rate, background_rate = u.Quantity(signal_rate, u.s**-1), u.Quantity(background_rate, u.s**-1)
    if signal_rate <= 0 or signal_rate < min_excess_over_background * background_rate:
        return np.inf * u.h
    alpha = 1 / n_off
    second = 1 * u.s
    n_off_counts = (background_rate * second * n_off).to_value(u.one)
    n_on = (signal_rate * second).to_value(u.one) + alpha * n_off_counts
    ts = WStatCountsStatistic(n_on=n_on, n_off=n_off_counts, alpha=alpha).ts
    time = max(sigma**2 / ts * second, min_signal_events / signal_rate)
    return time.to(u.h)


def time_to_detection_table(irf_file, sources, energies_min, **kwargs):
    """
    Time to detection of each source of ``sources`` (names of `sst1mpipe.performance.spectra`) for each
    minimum reconstructed energy of ``energies_min``: table (source, energy_min, signal and background
    rates, time); ``kwargs``: of `time_to_detection`
    """
    rows = []
    for source in sources:
        rates = detection_rates(irf_file, source_spectrum(source))
        for energy_min in u.Quantity(energies_min, u.TeV):
            signal, background = rates.above(energy_min)
            rows.append({"source": source, "energy_min": energy_min, "signal_rate": signal.to(u.min**-1),
                         "background_rate": background.to(u.min**-1),
                         "time": time_to_detection(signal, background, **kwargs)})
    return QTable(rows=rows)

