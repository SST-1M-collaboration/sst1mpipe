"""
Run a ctapipe IRF tool with the proton + helium spectrum of pyirf (DAMPE_P_He_SPECTRUM) instead of
IRFDOC_PROTON_SPECTRUM for the background (the protons of the simulations weighted to protons +
helium, as in the SST-1M Crab paper, arXiv:2506.01733). ctapipe has no option for it.

    python run_p_he.py ctapipe-optimize-event-selection <options>
    python run_p_he.py ctapipe-compute-irf <options>
"""
import sys

from ctapipe.irf import spectra
from pyirf.spectral import DAMPE_P_He_SPECTRUM

spectra.SPECTRA[spectra.Spectra.IRFDOC_PROTON_SPECTRUM] = DAMPE_P_He_SPECTRUM

TOOLS = {
    "ctapipe-optimize-event-selection": "ctapipe.tools.optimize_event_selection",
    "ctapipe-compute-irf": "ctapipe.tools.compute_irf",
}

if __name__ == "__main__":
    tool = sys.argv[1]
    module = __import__(TOOLS[tool], fromlist=["main"])
    sys.argv = [tool, *sys.argv[2:]]
    module.main()
