"""
IRFs with the ctapipe IRF tools as in the SST-1M Crab paper (arXiv:2506.01733):

- `use_spectra`: the spectra of the event weights of the ctapipe tools, the Crab Nebula of MAGIC
  (log-parabola of Aleksić et al. 2015, ``CRAB_MAGIC_JHEAP2015`` of pyirf, Sec. 5.1, 5.3 of the paper)
  in place of HEGRA, and the protons weighted to protons + helium (DAMPE, ``DAMPE_P_He_SPECTRUM``)
- `GlobalThetaPercentileCuts`: the cuts of the paper (Sec. 5.3), a gammaness cut keeping a fraction of
  the gammas in each reconstructed energy bin and one theta cut for all the energies, keeping a
  fraction of the gammas left by the gammaness cut
- `FixedThetaCuts`: the gammaness cut of `GlobalThetaPercentileCuts` and a fixed theta cut for all the
  energies (``theta_cut``, 0.2 deg: the cut of the analysis of the paper with sst1mpipe < 0.9,
  ``global_theta_cut``)
- `MinimumStatisticsSensitivity2dMaker`: the sensitivity of ctapipe (``Sensitivity2dMaker``) only in the
  energy bins with more than ``min_simulated_signal`` simulated gammas in the ON region and more than
  ``min_simulated_background`` simulated protons in the background region (10, as sst1mpipe < 0.9),
  NaN in the other ones

ctapipe has fixed names for its spectra (``Spectra``): ``CRAB_HEGRA`` and ``IRFDOC_PROTON_SPECTRUM`` of
the configurations are these spectra in the sst1mpipe tools (``sst1mpipe-optimize-event-selection``,
``sst1mpipe-compute-irf``).
"""
import operator

import astropy.units as u
import numpy as np
from astropy.io.fits import BinTableHDU
from astropy.table import QTable
from ctapipe.core.traits import AstroQuantity, Int
from ctapipe.irf import spectra
from ctapipe.irf.benchmarks import Sensitivity2dMaker
from ctapipe.irf.optimize import PercentileCuts
from pyirf.cuts import evaluate_binned_cut
from pyirf.spectral import CRAB_HEGRA, CRAB_MAGIC_JHEAP2015, DAMPE_P_He_SPECTRUM, IRFDOC_PROTON_SPECTRUM

__all__ = [
    "CRAB_SPECTRA",
    "BACKGROUND_SPECTRA",
    "FixedThetaCuts",
    "GlobalThetaPercentileCuts",
    "MinimumStatisticsSensitivity2dMaker",
    "use_spectra",
]

#: Crab Nebula spectrum of ``CRAB_HEGRA`` in the ctapipe IRF tools
CRAB_SPECTRA = {"MAGIC": CRAB_MAGIC_JHEAP2015, "HEGRA": CRAB_HEGRA}
#: proton spectrum of ``IRFDOC_PROTON_SPECTRUM`` in the ctapipe IRF tools: p_he (protons + helium), protons
BACKGROUND_SPECTRA = {"p_he": DAMPE_P_He_SPECTRUM, "protons": IRFDOC_PROTON_SPECTRUM}


def use_spectra(crab="MAGIC", background="p_he"):
    """
    Spectra of the event weights of the ctapipe IRF tools (``ctapipe.irf.spectra.SPECTRA``, global):
    ``CRAB_HEGRA`` -> the Crab spectrum ``crab`` (`CRAB_SPECTRA`), ``IRFDOC_PROTON_SPECTRUM`` -> the proton
    spectrum ``background`` (`BACKGROUND_SPECTRA`)
    """
    spectra.SPECTRA[spectra.Spectra.CRAB_HEGRA] = CRAB_SPECTRA[crab]
    spectra.SPECTRA[spectra.Spectra.IRFDOC_PROTON_SPECTRUM] = BACKGROUND_SPECTRA[background]


class GlobalThetaPercentileCuts(PercentileCuts):
    """
    Cuts of the SST-1M Crab paper: gammaness cut keeping ``GhPercentileCutCalculator.target_percentile`` %
    of the gammas in each reconstructed energy bin (as ``PercentileCuts``), and one theta cut for all the
    energies (the same value in every bin), keeping ``ThetaPercentileCutCalculator.target_percentile`` %
    of the gammas left by the gammaness cut in the energy range of the bins (events not weighted, as
    the cuts per bin), within ``theta_min_angle`` and ``theta_max_angle`` (``theta_fill_value`` with
    fewer than ``min_counts`` events)
    """

    def __call__(self, events, quality_query, clf_prefix):
        result = super().__call__(events, quality_query, clf_prefix)
        signal = events["signal"]
        gh_mask = evaluate_binned_cut(
            values=signal["gh_score"], bin_values=signal["reco_energy"], cut_table=result.gh_cuts, op=operator.ge,
        )
        energy = signal["reco_energy"]
        in_range = (energy >= self.reco_energy_bins[0]) & (energy < self.reco_energy_bins[-1])
        theta = signal["theta"][gh_mask & in_range].to(u.deg)
        calculator = self.theta_cut_calculator
        if len(theta) < calculator.min_counts:
            cut = calculator.theta_fill_value.to(u.deg)
        else:
            cut = np.percentile(theta, calculator.target_percentile)
            if calculator.theta_min_angle is not None:
                cut = max(cut, calculator.theta_min_angle.to(u.deg))
            if calculator.theta_max_angle is not None:
                cut = min(cut, calculator.theta_max_angle.to(u.deg))
        table = result.spatial_selection_table
        table["cut"] = u.Quantity(np.full(len(table), cut.to_value(u.deg)), u.deg).to(table["cut"].unit)
        self.log.info("theta cut for all the energies: %.3f deg (%d %%)", cut.to_value(u.deg),
                      calculator.target_percentile)
        return result


class FixedThetaCuts(PercentileCuts):
    """
    Gammaness cut keeping ``GhPercentileCutCalculator.target_percentile`` % of the gammas in each
    reconstructed energy bin (as ``PercentileCuts``), and a fixed theta cut ``theta_cut`` for all the
    energies (the analysis of the SST-1M Crab paper with sst1mpipe < 0.9: ``global_theta_cut`` 0.2 deg,
    about 80 % containment). Selected with ``EventSelectionOptimizer.optimization_algorithm:
    FixedThetaCuts`` and ``FixedThetaCuts.theta_cut``.
    """

    theta_cut = AstroQuantity(
        default_value=u.Quantity(0.2, u.deg),
        physical_type=u.physical.angle,
        help="Theta cut for all the energies",
    ).tag(config=True)

    def __call__(self, events, quality_query, clf_prefix):
        result = super().__call__(events, quality_query, clf_prefix)
        table = result.spatial_selection_table
        table["cut"] = u.Quantity(np.full(len(table), self.theta_cut.to_value(u.deg)), u.deg).to(table["cut"].unit)
        self.log.info("theta cut for all the energies: %.3f deg (fixed)", self.theta_cut.to_value(u.deg))
        return result
class MinimumStatisticsSensitivity2dMaker(Sensitivity2dMaker):
    """
    ``Sensitivity2dMaker`` of ctapipe (pyirf ``calculate_sensitivity``: 5 sigma, excess of at least 10
    weighted events and 5 % of the background) with the condition of sst1mpipe < 0.9 on the statistics of
    the simulations: the sensitivity of an energy bin (``RELATIVE_SENSITIVITY``, ``FLUX_SENSITIVITY``,
    ``ENERGY_FLUX_SENSITIVITY``, ``SIGNIFICANCE``) is NaN unless more than ``min_simulated_signal``
    simulated gammas are left in the ON region and more than ``min_simulated_background`` simulated
    protons in the background region (the offset bin, see ``fov_offset_*``), after the cuts. These
    numbers of simulated events (not weighted) are written in ``N_SIGNAL_SIMULATED`` and
    ``N_BACKGROUND_SIMULATED`` (``N_SIGNAL`` and ``N_BACKGROUND`` of ctapipe are scaled). Selected with
    ``IrfTool.sensitivity_maker_name: MinimumStatisticsSensitivity2dMaker``; the settings of the
    ``Sensitivity2dMaker`` section apply.
    """

    min_simulated_signal = Int(
        default_value=10, min=0,
        help="Sensitivity of an energy bin only with more than this number of simulated gammas in the ON region",
    ).tag(config=True)

    min_simulated_background = Int(
        default_value=10, min=0,
        help="Sensitivity of an energy bin only with more than this number of simulated protons in the background region",
    ).tag(config=True)

    def __call__(self, signal_events, background_events, spatial_selection_table, gamma_spectrum,
                 extname="SENSITIVITY"):
        from ctapipe.irf.benchmarks import _get_2d_result_table

        hdu = super().__call__(signal_events, background_events, spatial_selection_table, gamma_spectrum, extname)
        table = QTable.read(hdu)
        energy_bins = self.reco_energy_bins
        _, fov_bin_idx, shape = _get_2d_result_table(signal_events, energy_bins, self.fov_offset_bins)
        n_signal, n_background = np.zeros(shape, dtype=np.int64), np.zeros(shape, dtype=np.int64)
        offset = background_events["reco_source_fov_offset"]
        for i in range(len(self.fov_offset_bins) - 1):
            n_signal[i], _ = np.histogram(signal_events["reco_energy"][fov_bin_idx == i], energy_bins)
            in_ring = (offset >= self.fov_offset_bins[i]) & (offset < self.fov_offset_bins[i + 1])
            n_background[i], _ = np.histogram(background_events["reco_energy"][in_ring], energy_bins)
        invalid = (n_signal <= self.min_simulated_signal) | (n_background <= self.min_simulated_background)
        for name in ("SIGNIFICANCE", "RELATIVE_SENSITIVITY", "FLUX_SENSITIVITY", "ENERGY_FLUX_SENSITIVITY"):
            values = table[name].copy()
            values[0][invalid] = np.nan
            table[name] = values
        table["N_SIGNAL_SIMULATED"] = n_signal[np.newaxis]
        table["N_BACKGROUND_SIMULATED"] = n_background[np.newaxis]
        self.log.info(
            "sensitivity: %d / %d energy bins with more than %d simulated gammas and %d simulated protons",
            np.count_nonzero(~invalid), invalid.size, self.min_simulated_signal, self.min_simulated_background,
        )
        header = hdu.header.copy()
        header["MINSIMS"] = (self.min_simulated_signal, "Simulated gammas per bin: more than")
        header["MINSIMB"] = (self.min_simulated_background, "Simulated protons per bin: more than")
        return BinTableHDU(table, header=header, name=extname)
