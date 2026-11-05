import numpy as np
from ctapipe.instrument.subarray import SubarrayDescription
from sst1mpipe.instrument.camera import Camera
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
# They are built from the camera configuration file read by CTS (camera_config.cfg):
#
# - ``trigger_input_traces`` (FADC output of the patches, 192 channels per sector of which the
#   first 144 are used): position 144 * (micro_crate - 1) + 4 * rank of module_fw_id in the
#   sector (36 modules) + patch_in_mod_fw
# - ``trigger_output_patch7``, ``trigger_output_patch19`` and ``trigger_output_muon`` (one bit per
#   cluster, 18 bytes per sector): position 144 * (micro_crate - 1) + patch_trig_out

N_SECTORS = 3
PATCHES_PER_SECTOR = 144
N_TRACES_PER_SECTOR = 192  # trigger input traces sent per sector, the first PATCHES_PER_SECTOR are the patches
N_PATCHES = N_SECTORS * PATCHES_PER_SECTOR  # 432

# columns of the camera configuration file (one row per pixel) used for the trigger patches
CAMERA_CONFIG_TRIGGER_COLUMNS = {
    "micro_crate": 3,
    "patch_sw_id": 11,
    "patch_in_mod_fw": 13,
    "module_fw_id": 14,
    "patch_trig_out": 15,
}


def read_trigger_patches(config_file=CAMERA.config_file):
    """
    Readout information of the trigger patches, from the camera configuration file (one row
    per pixel, the same values for the 3 pixels of a patch).

    Parameters
    ----------
    config_file: path
        camera_config.cfg, by default the one of the camera (`CAMERA`)

    Returns
    -------
    dict:
        column name (`CAMERA_CONFIG_TRIGGER_COLUMNS`) -> np.ndarray (n_patches),
        indexed by patch_sw_id
    """
    columns = list(CAMERA_CONFIG_TRIGGER_COLUMNS)
    rows = np.loadtxt(config_file, comments="#", usecols=list(CAMERA_CONFIG_TRIGGER_COLUMNS.values()), dtype=int)
    rows = rows[np.argsort(rows[:, columns.index("patch_sw_id")], kind="stable")]
    patch_ids = rows[:, columns.index("patch_sw_id")]
    patches, first = np.unique(patch_ids, return_index=True)
    if not np.array_equal(patches, np.arange(len(patches))):
        raise ValueError(f"The patch_sw_id of {config_file} are not 0 to n_patches - 1")
    patch_rows = rows[first]
    if not np.array_equal(rows, patch_rows[patch_ids]):
        raise ValueError(f"The pixels of a patch have different readout values in {config_file}")
    return {name: patch_rows[:, i] for i, name in enumerate(columns)}


def build_patch_ids(patches=None):
    """
    patch_sw_id of each position of the trigger traces for the sectors (1, 2, 3), see above

    Parameters
    ----------
    patches: dict or None
        Readout information of the patches (`read_trigger_patches`), by default of the camera

    Returns
    -------
    dict: ``PATCH_IDS``, kind of trigger traces ("input", "output") -> list (n_patches)
    """
    patches = read_trigger_patches() if patches is None else patches
    sector_start = PATCHES_PER_SECTOR * (patches["micro_crate"] - 1)
    patches_per_module = patches["patch_in_mod_fw"].max() + 1
    module_rank = np.empty_like(patches["module_fw_id"])
    for sector in np.unique(patches["micro_crate"]):
        in_sector = patches["micro_crate"] == sector
        module_rank[in_sector] = np.unique(patches["module_fw_id"][in_sector], return_inverse=True)[1]

    def patch_ids(position):
        ids = np.full(N_PATCHES, -1)
        ids[position] = np.arange(len(position))
        if np.any(ids < 0):
            raise ValueError("Two trigger patches at the same position of the trigger traces")
        return ids.tolist()

    input_position = sector_start + patches_per_module * module_rank + patches["patch_in_mod_fw"]
    output_position = sector_start + patches["patch_trig_out"]
    return {"input": patch_ids(input_position), "output": patch_ids(output_position)}


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
