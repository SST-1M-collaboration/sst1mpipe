"""
Machine learning reconstruction of the SST-1M events with the ctapipe reconstructors
(``ctapipe.reco``: energy regressor, particle classifier, disp reconstructor), trained and applied
by the sst1mpipe tools (``sst1mpipe.scripts``): the stereo combiners of the direction from the disp
norms and a choice of the disp signs (`StereoDispSignClosestCombiner`: SST-1M Crab paper,
`StereoDispSignCrossingCombiner`: crossing of the major axes) and the plots of the trained models
(`sst1mpipe.ml.plots`). Importing this package makes the combiners available to the ctapipe
reconstructors.
"""
from .stereo_combination import (
    StereoDispSignClosestCombiner,
    StereoDispSignCombiner,
    StereoDispSignCrossingCombiner,
    axes_crossing_disp_signs,
    disp_sign_combination,
)

__all__ = [
    "StereoDispSignCombiner",
    "StereoDispSignClosestCombiner",
    "StereoDispSignCrossingCombiner",
    "disp_sign_combination",
    "axes_crossing_disp_signs",
]
