# Licensed under a 3-clause BSD style license - see LICENSE
"""
Performance of the simulations: weights of the simulated events (`weights`), spectra of sources
(`spectra`), figures of the IRFs and benchmarks of ctapipe and of reference instruments (`irf_plots`)
"""

from .weights import get_mc_info, get_weights

__all__ = ["get_mc_info", "get_weights"]
