import astropy.units as u
import numpy as np
from ctapipe.instrument import SubarrayDescription
from cts_core.camera import Camera

from sst1mpipe.constants import (N_PIXELS, CAMERA, SUBARRAY_DESCRIPTION, N_CHANNELS,
                                 PATCH_ID_INPUT, PATCH_ID_INPUT_SORT_IDS, PATCH_ID_OUTPUT, PATCH_ID_OUTPUT_SORT_IDS)

def test_constants():

    assert isinstance(N_PIXELS, int)
    assert isinstance(CAMERA, Camera)
    assert isinstance(SUBARRAY_DESCRIPTION, SubarrayDescription)
    assert isinstance(N_CHANNELS, int)
    assert isinstance(PATCH_ID_INPUT,list)
    assert isinstance(PATCH_ID_OUTPUT, list)
    assert isinstance(PATCH_ID_OUTPUT_SORT_IDS, np.ndarray)
    assert isinstance(PATCH_ID_INPUT_SORT_IDS, np.ndarray)


def test_subarray_reference_location():
    # the reference location of the subarray is the SST-1M site (Ondrejov): it is used by ctapipe
    # (GroundFrame, AltAz of the pointing...), there is no other reference location in sst1mpipe
    location = SUBARRAY_DESCRIPTION.reference_location.geodetic
    assert np.isclose(location.lon.to_value(u.deg), 14.7825, atol=1e-3)
    assert np.isclose(location.lat.to_value(u.deg), 49.9124, atol=1e-3)
    assert np.isclose(location.height.to_value(u.m), 511.2, atol=1)
