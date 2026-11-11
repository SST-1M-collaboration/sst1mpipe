import numpy as np
from astropy import units as u
from ctapipe.coordinates import EngineeringCameraFrame
from ctapipe.instrument.camera import CameraGeometry
from scipy.sparse import identity

# from ctapipe.instrument.camera import find_neighbor_pixels
from sst1m_cts_core import camera as cam

from sst1mpipe.resources import CAMERA_CONFIG_FILE


def find_pixel_positions(
        camera_config_file,
        source_x=0. * u.mm,
        source_y=0. * u.mm
):
    camera = cam.Camera(_config_file=camera_config_file)

    x, y = [], []

    for pixel in camera.Pixels:
        x.append(pixel.center[0])
        y.append(pixel.center[1])

    r = np.array([x, y]) * u.mm
    r[0] = r[0] - source_x
    r[1] = r[1] - source_y

    return r


def generate_geometry_from_camera(
        camera,
        source_x=0. * u.mm,
        source_y=0. * u.mm,
        name = 'SST-1M'
    ):
    """
    Generate the SST-1M geometry from the CTS configuration.

    The pixel positions of the CTS configuration are in the ctapipe `EngineeringCameraFrame`
    (camera seen from the front, x = -y and y = -x of the `CameraFrame`): the geometry has this
    frame, so that ctapipe draws it in this orientation and transforms it to the `CameraFrame` of the
    subarray of the event sources and of the simulations (``geometry.transform_to(CameraFrame())``).

    :param camera: a CTS instance
    :param source_x:  x coord
    :param source_y:  y coord
    :return: the geometry for visualisation and the list of "good" pixels
    """
    pix_x = []
    pix_y = []
    pix_area = []
    pix_id = []

    for pix in camera.Pixels:
        pix_x.append(pix.center[0])
        pix_y.append(pix.center[1])
        pix_id.append(pix.ID)
        pix_area.append(pix.area)

    # neighbors_pix = _find_neighbor_pixels(pix_x, pix_y, 30.)
    pix_x = np.array(pix_x) * u.mm - source_x
    pix_y = np.array(pix_y) * u.mm - source_y
    pix_area = np.array(pix_area) * u.mm**2

    geom = CameraGeometry(
        name,
        pix_id,
        pix_x,
        pix_y,
        pix_area,
        pix_type='hexagonal',
        frame=EngineeringCameraFrame(),
        # neighbors=neighbors_pix
    )

    return geom


#: type of the DigiCam matrices (patches, clusters): their products stay exact (at most
#: 19 patches x 255 = 4845) and scipy keeps int16 for the products with the uint8 trigger traces,
#: about 4 times faster than int64 (bool or uint8 matrices would give uint8 products, which overflow)
MATRIX_DTYPE = np.int16


def compute_patch_matrix(camera):
    n_pixels = len(camera.Pixels)
    n_patches = len(camera.Patches)
    patch_matrix = np.zeros((n_patches, n_pixels), dtype=MATRIX_DTYPE)

    for patch in camera.Patches:

        for pixel in patch.pixelsID:
            patch_matrix[patch.ID, pixel] = 1

    return patch_matrix


def compute_cluster_matrix_7(camera):
    n_clusters = len(camera.Clusters_7)
    cluster_matrix = np.zeros((n_clusters, n_clusters), dtype=MATRIX_DTYPE)

    for cluster in camera.Clusters_7:

        for patch in cluster.patchesID:
            cluster_matrix[cluster.ID, patch] = 1

    return cluster_matrix


def compute_neighbor_patch_matrix(cluster_7_matrix):
    """
    (n_patches, n_patches) sparse matrix, 1 where two patches are neighbours: the patch7 clusters
    without their centre
    """
    n_patches = cluster_7_matrix.shape[0]
    return (cluster_7_matrix - identity(n_patches, dtype=MATRIX_DTYPE, format="csr")).tocsr()


def compute_cluster_matrix_19(camera):
    n_clusters = len(camera.Clusters_19)
    cluster_matrix = np.zeros((n_clusters, n_clusters), dtype=MATRIX_DTYPE)

    for cluster in camera.Clusters_19:

        for patch in cluster.patchesID:
            cluster_matrix[cluster.ID, patch] = 1

    return cluster_matrix


# columns of the camera configuration file (one row per pixel) used for the trigger patches
CAMERA_CONFIG_TRIGGER_COLUMNS = {
    "micro_crate": 3,
    "patch_sw_id": 11,
    "patch_in_mod_fw": 13,
    "module_fw_id": 14,
    "patch_trig_out": 15,
}


def read_trigger_patches(config_file=CAMERA_CONFIG_FILE):
    """
    Readout information of the trigger patches, from the camera configuration file (one row
    per pixel, the same values for the 3 pixels of a patch).

    Parameters
    ----------
    config_file: path
        camera_config.cfg, by default the one of sst1mpipe (`CAMERA_CONFIG_FILE`)

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
    patch_sw_id of each position of the trigger traces of the camera events when their blocks are
    the sectors (1, 2, 3).

    The camera is read by sectors (micro_crate) of patches (3 pixels, patch_sw_id), and the trigger
    traces are made of one block of patches per sector:

    - ``trigger_input_traces`` (FADC output of the patches): position
      patches_per_sector * (micro_crate - 1) + patches_per_module * rank of module_fw_id in the
      sector + patch_in_mod_fw
    - ``trigger_output_patch7``, ``trigger_output_patch19`` and ``trigger_output_muon`` (one bit
      per cluster): position patches_per_sector * (micro_crate - 1) + patch_trig_out

    Parameters
    ----------
    patches: dict or None
        Readout information of the patches (`read_trigger_patches`), by default of the camera

    Returns
    -------
    dict: ``PATCH_IDS``, kind of trigger traces ("input", "output") -> list (n_patches)
    """
    patches = read_trigger_patches() if patches is None else patches
    n_patches = len(patches["micro_crate"])
    patches_per_sector = n_patches // len(np.unique(patches["micro_crate"]))
    sector_start = patches_per_sector * (patches["micro_crate"] - 1)
    patches_per_module = patches["patch_in_mod_fw"].max() + 1
    module_rank = np.empty_like(patches["module_fw_id"])
    for sector in np.unique(patches["micro_crate"]):
        in_sector = patches["micro_crate"] == sector
        module_rank[in_sector] = np.unique(patches["module_fw_id"][in_sector], return_inverse=True)[1]

    def patch_ids(position):
        ids = np.full(n_patches, -1)
        ids[position] = np.arange(len(position))
        if np.any(ids < 0):
            raise ValueError("Two trigger patches at the same position of the trigger traces")
        return ids.tolist()

    input_position = sector_start + patches_per_module * module_rank + patches["patch_in_mod_fw"]
    output_position = sector_start + patches["patch_trig_out"]
    return {"input": patch_ids(input_position), "output": patch_ids(output_position)}
