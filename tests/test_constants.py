import astropy.units as u
import numpy as np
from ctapipe.instrument import SubarrayDescription
from cts_core.camera import Camera

from sst1mpipe.constants import (N_PIXELS, CAMERA, SUBARRAY_DESCRIPTION, N_CHANNELS, N_PATCHES, PATCHES_PER_SECTOR,
                                 PATCH_IDS, PATCH_SORT_IDS,
                                 TRIGGER_SECTOR_ORDER, patch_sort_ids)
from sst1mpipe.instrument.geometry import build_patch_ids, read_trigger_patches

def test_constants():

    assert isinstance(N_PIXELS, int)
    assert isinstance(CAMERA, Camera)
    assert isinstance(SUBARRAY_DESCRIPTION, SubarrayDescription)
    assert isinstance(N_CHANNELS, int)
    assert set(PATCH_IDS) == set(PATCH_SORT_IDS) == {"input", "output"}
    assert all(isinstance(patch_ids, list) for patch_ids in PATCH_IDS.values())
    assert all(isinstance(sort_ids, dict) for sort_ids in PATCH_SORT_IDS.values())


def test_subarray_reference_location():
    # the reference location of the subarray is the SST-1M site (Ondrejov): it is used by ctapipe
    # (GroundFrame, AltAz of the pointing...), there is no other reference location in sst1mpipe
    location = SUBARRAY_DESCRIPTION.reference_location.geodetic
    assert np.isclose(location.lon.to_value(u.deg), 14.7825, atol=1e-3)
    assert np.isclose(location.lat.to_value(u.deg), 49.9124, atol=1e-3)
    assert np.isclose(location.height.to_value(u.m), 511.2, atol=1)


def test_patch_ids_of_the_camera():
    assert PATCH_IDS == build_patch_ids()


def test_patch_sort_ids_of_the_telescopes():
    patches = read_trigger_patches()
    for kind, patch_ids in PATCH_IDS.items():
        sort_ids = PATCH_SORT_IDS[kind]
        # sectors (1, 2, 3), as tel 21: the sorted patch ids
        default = patch_sort_ids(patch_ids, (1, 2, 3))
        assert set(sort_ids) == set(TRIGGER_SECTOR_ORDER) == {21, 22}
        assert np.array_equal(np.asarray(patch_ids)[default], np.arange(N_PATCHES))
        assert np.array_equal(sort_ids[21], default)
        for tel_id, sector_order in TRIGGER_SECTOR_ORDER.items():
            # position of each patch in the traces of the telescope: in the block of its sector,
            # at the same place inside the block for all the telescopes
            position = sort_ids[tel_id]
            assert sorted(position) == list(range(N_PATCHES))
            block = np.array([sector_order.index(sector) for sector in patches["micro_crate"]])
            assert np.array_equal(position // PATCHES_PER_SECTOR, block)
            assert np.array_equal(position % PATCHES_PER_SECTOR, default % PATCHES_PER_SECTOR)
    # telescope 2: sectors shifted by one
    assert not np.array_equal(PATCH_SORT_IDS["input"][22], PATCH_SORT_IDS["input"][21])
    assert np.array_equal(patch_sort_ids(PATCH_IDS["input"], (2, 3, 1)), PATCH_SORT_IDS["input"][22])
