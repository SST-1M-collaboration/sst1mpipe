import numpy as np
from ctapipe.instrument.subarray import SubarrayDescription
from sst1mpipe.instrument.camera import Camera
from sst1mpipe.instrument.geometry import build_patch_ids
from sst1mpipe.resources import SUBARRAY_FILE

SUBARRAY_DESCRIPTION = SubarrayDescription.from_hdf(SUBARRAY_FILE)

CAMERA = Camera()
GEOMETRY = CAMERA.geometry
N_PIXELS = GEOMETRY.n_pixels
N_CHANNELS = 1

# ---------------------------------------------------------------------------
# Order of the trigger traces of the DigiCam camera events
# ---------------------------------------------------------------------------
# The camera is read by three sectors (micro_crate in camera_config.cfg) of 144 trigger patches
# (patch of 3 pixels, patch_sw_id). The trigger traces of the camera events are made of three
# blocks of 144 patches, one per sector: ``PATCH_IDS["input"]`` and ``PATCH_IDS["output"]`` give
# the patch_sw_id of each position of the traces when the blocks are the sectors (1, 2, 3).
# They are built from the camera configuration file read by CTS (camera_config.cfg), see
# `sst1mpipe.instrument.geometry.build_patch_ids`.

N_SECTORS = 3
PATCHES_PER_SECTOR = 144
N_TRACES_PER_SECTOR = 192  # trigger input traces sent per sector, the first PATCHES_PER_SECTOR are the patches
N_PATCHES = N_SECTORS * PATCHES_PER_SECTOR  # 432

# kind of trigger traces ("input": trigger_input_traces, "output": the trigger outputs) ->
# patch_sw_id of each position of the traces for the sectors (1, 2, 3)
PATCH_IDS = build_patch_ids()

# Order of the sectors in the blocks of the trigger traces of each telescope. Measured by comparing
# trigger_input_traces with the patches rebuilt from the pixel waveforms (correlation 0.997 for
# tel 21, 0.999 for tel 22), and confirmed on trigger_output_patch7 (identical to the emulated
# patch7 output, sample by sample, on both telescopes).
TRIGGER_SECTOR_ORDER = {21: (1, 2, 3), 22: (2, 3, 1)}


def patch_sort_ids(patch_ids, sector_order):
    """
    Indices putting trigger traces in patch_sw_id order: ``traces[sort_ids][p]`` is the patch ``p``.

    Parameters
    ----------
    patch_ids: list
        patch_sw_id of each position of the traces for the sectors (1, 2, 3), see ``PATCH_IDS``
    sector_order: tuple
        Sector (micro_crate) of each block of 144 patches of the traces of the telescope
    """
    blocks = np.reshape(patch_ids, (N_SECTORS, PATCHES_PER_SECTOR))
    return np.argsort(blocks[np.subtract(sector_order, 1)].ravel())


# kind of trigger traces -> tel_id -> indices putting the traces of the telescope in patch_sw_id order
PATCH_SORT_IDS = {
    kind: {tel_id: patch_sort_ids(patch_ids, order) for tel_id, order in TRIGGER_SECTOR_ORDER.items()}
    for kind, patch_ids in PATCH_IDS.items()
}
