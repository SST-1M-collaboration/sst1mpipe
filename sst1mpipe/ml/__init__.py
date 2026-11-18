"""
Machine learning reconstruction of the SST-1M events with the ctapipe reconstructors
(``ctapipe.reco``: energy regressor, particle classifier, disp reconstructor), trained and applied
by the sst1mpipe tools (``sst1mpipe.scripts``): `StereoMARSCombiner` (stereo direction of the MARS
analysis) and the plots of the trained models. Importing this package makes `StereoMARSCombiner`
available to the ctapipe reconstructors.
"""
from .stereo_combination import StereoMARSCombiner, mars_directions

__all__ = ["StereoMARSCombiner", "mars_directions"]
