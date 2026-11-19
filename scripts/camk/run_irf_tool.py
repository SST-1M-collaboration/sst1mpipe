"""
Run a ctapipe IRF tool (ctapipe-optimize-event-selection, ctapipe-compute-irf) with the background
spectrum IRF_BACKGROUND (environment variable):

- ``p_he`` (default): the protons of the simulations weighted to protons + helium
  (``DAMPE_P_He_SPECTRUM`` of pyirf, as the SST-1M Crab paper, arXiv:2506.01733), in place of
  ``IRFDOC_PROTON_SPECTRUM`` (ctapipe has no option for it)
- ``protons``: ``IRFDOC_PROTON_SPECTRUM`` of ctapipe

    python run_irf_tool.py ctapipe-compute-irf <options>
"""
import os
import sys

from ctapipe.irf import spectra
from pyirf.spectral import DAMPE_P_He_SPECTRUM

TOOLS = {
    "ctapipe-optimize-event-selection": "ctapipe.tools.optimize_event_selection",
    "ctapipe-compute-irf": "ctapipe.tools.compute_irf",
}

if __name__ == "__main__":
    background = os.environ.get("IRF_BACKGROUND", "p_he")
    if background == "p_he":
        spectra.SPECTRA[spectra.Spectra.IRFDOC_PROTON_SPECTRUM] = DAMPE_P_He_SPECTRUM
    elif background != "protons":
        sys.exit(f"IRF_BACKGROUND: p_he or protons, not {background}")
    tool = sys.argv[1]
    module = __import__(TOOLS[tool], fromlist=["main"])
    sys.argv = [tool, *sys.argv[2:]]
    module.main()
