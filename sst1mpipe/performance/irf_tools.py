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

ctapipe has fixed names for its spectra (``Spectra``): ``CRAB_HEGRA`` and ``IRFDOC_PROTON_SPECTRUM`` of
the configurations are these spectra in the sst1mpipe tools (``sst1mpipe-optimize-event-selection``,
``sst1mpipe-compute-irf``).
"""
import operator

import astropy.units as u
import numpy as np
from ctapipe.irf import spectra
from ctapipe.core.traits import AstroQuantity
from ctapipe.irf.optimize import PercentileCuts
from pyirf.cuts import evaluate_binned_cut
from pyirf.spectral import CRAB_HEGRA, CRAB_MAGIC_JHEAP2015, DAMPE_P_He_SPECTRUM, IRFDOC_PROTON_SPECTRUM

__all__ = ["CRAB_SPECTRA", "BACKGROUND_SPECTRA", "FixedThetaCuts", "GlobalThetaPercentileCuts", "use_spectra"]

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
