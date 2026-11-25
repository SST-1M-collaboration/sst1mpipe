"""
Machine learning reconstruction of the SST-1M events with the ctapipe reconstructors
(``ctapipe.reco``: energy regressor, particle classifier, disp reconstructor), trained and applied
by the sst1mpipe tools (``sst1mpipe.scripts``): `StereoDISPSignCombiner` (stereo direction of the
SST-1M Crab paper) and the plots of the trained models (`sst1mpipe.ml.plots`). Importing this package
makes `StereoDISPSignCombiner` available to the ctapipe reconstructors.
"""
from .stereo_combination import StereoDISPSignCombiner, axes_crossing_disp_signs, disp_sign_combination

__all__ = ["StereoDISPSignCombiner", "disp_sign_combination", "axes_crossing_disp_signs"]
