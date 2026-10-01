import numpy as np
from astropy.coordinates import EarthLocation
from ctapipe.instrument import SubarrayDescription
from cts_core.camera import Camera

from sst1mpipe.constants import (N_PIXELS, CAMERA, REFERENCE_LOCATION, SUBARRAY_DESCRIPTION, N_CHANNELS,
                                 PATCH_ID_INPUT, PATCH_ID_INPUT_SORT_IDS, PATCH_ID_OUTPUT, PATCH_ID_OUTPUT_SORT_IDS)

def test_constants():

    assert isinstance(N_PIXELS, int)
    assert isinstance(CAMERA, Camera)
    assert isinstance(REFERENCE_LOCATION, EarthLocation)
    assert isinstance(SUBARRAY_DESCRIPTION, SubarrayDescription)
    assert isinstance(N_CHANNELS, int)
    assert isinstance(PATCH_ID_INPUT,list)
    assert isinstance(PATCH_ID_OUTPUT, list)
    assert isinstance(PATCH_ID_OUTPUT_SORT_IDS, np.ndarray)
    assert isinstance(PATCH_ID_INPUT_SORT_IDS, np.ndarray)

