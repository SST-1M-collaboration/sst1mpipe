import astropy.units as u
import numpy as np
from ctapipe.instrument import SubarrayDescription
from cts_core.camera import Camera

from sst1mpipe.constants import (N_PIXELS, CAMERA, SUBARRAY_DESCRIPTION, N_CHANNELS, N_PATCHES, PATCHES_PER_SECTOR,
                                 PATCH_ID_INPUT, PATCH_ID_INPUT_SORT_IDS, PATCH_ID_OUTPUT, PATCH_ID_OUTPUT_SORT_IDS,
                                 DEFAULT_PATCH_ID_INPUT_SORT_IDS, DEFAULT_PATCH_ID_OUTPUT_SORT_IDS,
                                 TRIGGER_SECTOR_ORDER, patch_sort_ids)
from sst1mpipe.resources import CAMERA_CONFIG_FILE

def test_constants():

    assert isinstance(N_PIXELS, int)
    assert isinstance(CAMERA, Camera)
    assert isinstance(SUBARRAY_DESCRIPTION, SubarrayDescription)
    assert isinstance(N_CHANNELS, int)
    assert isinstance(PATCH_ID_INPUT,list)
    assert isinstance(PATCH_ID_OUTPUT, list)
    assert isinstance(PATCH_ID_OUTPUT_SORT_IDS, dict)
    assert isinstance(PATCH_ID_INPUT_SORT_IDS, dict)


def test_subarray_reference_location():
    # the reference location of the subarray is the SST-1M site (Ondrejov): it is used by ctapipe
    # (GroundFrame, AltAz of the pointing...), there is no other reference location in sst1mpipe
    location = SUBARRAY_DESCRIPTION.reference_location.geodetic
    assert np.isclose(location.lon.to_value(u.deg), 14.7825, atol=1e-3)
    assert np.isclose(location.lat.to_value(u.deg), 49.9124, atol=1e-3)
    assert np.isclose(location.height.to_value(u.m), 511.2, atol=1)


def camera_config_patches():
    """camera_config.cfg columns of each patch (patch_sw_id -> dict), the same for its 3 pixels"""
    columns = {"micro_crate": 3, "patch_sw_id": 11, "patch_in_mod_fw": 13, "module_fw_id": 14, "patch_trig_out": 15}
    patches = {}
    with open(CAMERA_CONFIG_FILE) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            row = line.split()
            values = {name: int(row[column]) for name, column in columns.items()}
            assert patches.setdefault(values["patch_sw_id"], values) == values
    return patches


def test_patch_ids_from_camera_config():
    # PATCH_ID_INPUT and PATCH_ID_OUTPUT: patch_sw_id of each position of the trigger traces,
    # sectors (micro_crate) in the order (1, 2, 3), see sst1mpipe/constants.py
    patches = camera_config_patches()
    assert N_PATCHES == len(patches) == 432

    def patch_ids(position):
        ids = np.full(N_PATCHES, -1)
        for patch_id, patch in patches.items():
            ids[position(patch)] = patch_id
        assert sorted(ids) == list(range(N_PATCHES))
        return ids.tolist()

    module_rank = {module: i for i, module in enumerate(sorted({p["module_fw_id"] for p in patches.values()}))}
    sector = lambda patch: PATCHES_PER_SECTOR * (patch["micro_crate"] - 1)  # noqa: E731
    assert PATCH_ID_INPUT == patch_ids(
        lambda patch: sector(patch) + 4 * module_rank[patch["module_fw_id"]] + patch["patch_in_mod_fw"]
    )
    assert PATCH_ID_OUTPUT == patch_ids(lambda patch: sector(patch) + patch["patch_trig_out"])


def test_patch_sort_ids_of_the_telescopes():
    patches = camera_config_patches()
    for patch_ids, sort_ids, default in [
        (PATCH_ID_INPUT, PATCH_ID_INPUT_SORT_IDS, DEFAULT_PATCH_ID_INPUT_SORT_IDS),
        (PATCH_ID_OUTPUT, PATCH_ID_OUTPUT_SORT_IDS, DEFAULT_PATCH_ID_OUTPUT_SORT_IDS),
    ]:
        assert set(sort_ids) == set(TRIGGER_SECTOR_ORDER) == {21, 22}
        # sectors (1, 2, 3): the sorted patch ids
        assert np.array_equal(np.asarray(patch_ids)[default], np.arange(N_PATCHES))
        assert np.array_equal(sort_ids[21], default)
        for tel_id, sector_order in TRIGGER_SECTOR_ORDER.items():
            # position of each patch in the traces of the telescope: the block of its sector
            position = sort_ids[tel_id]
            assert sorted(position) == list(range(N_PATCHES))
            for patch_id, patch in patches.items():
                block = sector_order.index(patch["micro_crate"])
                assert position[patch_id] // PATCHES_PER_SECTOR == block
                # same place inside the block for all the telescopes
                assert position[patch_id] % PATCHES_PER_SECTOR == default[patch_id] % PATCHES_PER_SECTOR
    # telescope 2: sectors shifted by one
    assert not np.array_equal(PATCH_ID_INPUT_SORT_IDS[22], DEFAULT_PATCH_ID_INPUT_SORT_IDS)
    assert np.array_equal(patch_sort_ids(PATCH_ID_INPUT, (2, 3, 1)), PATCH_ID_INPUT_SORT_IDS[22])
