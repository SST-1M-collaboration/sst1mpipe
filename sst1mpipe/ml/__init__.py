"""
Machine learning reconstruction of the SST-1M events with the ctapipe reconstructors
(``ctapipe.reco``: energy regressor, particle classifier, disp reconstructor), trained and applied
by the sst1mpipe tools built on the ctapipe tools:

- ``sst1mpipe-train-energy-regressor``, ``sst1mpipe-train-particle-classifier``,
  ``sst1mpipe-train-disp-reconstructor``: the ctapipe training tools, with plots of the feature
  importances and of the performance of the models (cross validation);
- ``sst1mpipe-apply-models``: ``ctapipe-apply-models``, with the models chosen in a production of
  models per zenith, azimuth and NSB level;
- `StereoMARSCombiner`: the stereo direction of the MARS analysis, selected in the configuration
  of the disp reconstructor.

Importing this package makes `StereoMARSCombiner` available to the ctapipe reconstructors.
"""
from .stereo_combination import StereoMARSCombiner, mars_directions

__all__ = ["StereoMARSCombiner", "mars_directions"]
