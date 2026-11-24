"""
``sst1mpipe-optimize-event-selection`` and ``sst1mpipe-compute-irf``: ``ctapipe-optimize-event-selection``
and ``ctapipe-compute-irf`` with the spectra and cuts of the SST-1M Crab paper (`sst1mpipe.performance.irf_tools`):

- ``--crab MAGIC`` (default) or ``HEGRA``: Crab Nebula spectrum of ``CRAB_HEGRA`` in the configurations
- ``--background p_he`` (default, protons weighted to protons + helium) or ``protons``: proton spectrum of
  ``IRFDOC_PROTON_SPECTRUM``
- the cut algorithm ``GlobalThetaPercentileCuts`` (``sst1mpipe_irf_paper_cuts.yaml``)

The other options and configurations are those of the ctapipe tools::

    sst1mpipe-optimize-event-selection --config sst1mpipe_optimize_cuts.yaml --config sst1mpipe_irf_paper_cuts.yaml \\
        -g gamma_point.dl2.h5 -p proton_test.dl2.h5 -o cuts.fits
    sst1mpipe-compute-irf --config sst1mpipe_compute_irf.yaml --cuts cuts.fits --gamma-file gamma_point.dl2.h5 \\
        --proton-file proton_test.dl2.h5 --do-background --spatial-selection-applied -o irf.fits.gz
"""
from ctapipe.core.traits import CaselessStrEnum
from ctapipe.tools.compute_irf import IrfTool
from ctapipe.tools.optimize_event_selection import EventSelectionOptimizer

from sst1mpipe.performance.irf_tools import BACKGROUND_SPECTRA, CRAB_SPECTRA, GlobalThetaPercentileCuts, use_spectra

__all__ = ["SST1MEventSelectionOptimizer", "SST1MIrfTool", "GlobalThetaPercentileCuts"]

def _crab_trait():
    return CaselessStrEnum(
        list(CRAB_SPECTRA), default_value="MAGIC",
        help="Crab Nebula spectrum of CRAB_HEGRA: MAGIC (log-parabola, SST-1M Crab paper) or HEGRA",
    ).tag(config=True)


def _background_trait():
    return CaselessStrEnum(
        list(BACKGROUND_SPECTRA), default_value="p_he",
        help="Proton spectrum of IRFDOC_PROTON_SPECTRUM: p_he (protons + helium of DAMPE) or protons",
    ).tag(config=True)


class SST1MEventSelectionOptimizer(EventSelectionOptimizer):
    """ctapipe-optimize-event-selection with the spectra and cuts of the SST-1M Crab paper"""

    name = "sst1mpipe-optimize-event-selection"
    crab = _crab_trait()
    background = _background_trait()
    aliases = {**EventSelectionOptimizer.aliases, "crab": "SST1MEventSelectionOptimizer.crab",
               "background": "SST1MEventSelectionOptimizer.background"}

    def setup(self):
        use_spectra(crab=self.crab.upper(), background=self.background.lower())
        super().setup()


class SST1MIrfTool(IrfTool):
    """ctapipe-compute-irf with the spectra of the SST-1M Crab paper"""

    name = "sst1mpipe-compute-irf"
    crab = _crab_trait()
    background = _background_trait()
    aliases = {**IrfTool.aliases, "crab": "SST1MIrfTool.crab", "background": "SST1MIrfTool.background"}

    def setup(self):
        use_spectra(crab=self.crab.upper(), background=self.background.lower())
        super().setup()


def optimize_main():
    SST1MEventSelectionOptimizer().run()


def compute_irf_main():
    SST1MIrfTool().run()
