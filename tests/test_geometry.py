import numpy as np
import pytest

from sst1mpipe.constants import N_PATCHES, PATCHES_PER_SECTOR
from sst1mpipe.instrument.camera import Camera
from sst1mpipe.instrument.geometry import CAMERA_CONFIG_TRIGGER_COLUMNS, build_patch_ids, read_trigger_patches
from sst1mpipe.resources import CAMERA_CONFIG_FILE

CAMERA = Camera()
GEOMETRY = CAMERA.geometry

def test_geometry_pixel_has_units():

    assert GEOMETRY.pix_x.unit.physical_type == "length"
    assert GEOMETRY.pix_y.unit.physical_type == "length"
    assert GEOMETRY.pix_area.unit.physical_type == "area"

def test_get_neighbors_sparse_matrix():

    assert isinstance(GEOMETRY.neighbor_matrix, np.ndarray)


# sst1mpipe.constants.PATCH_IDS (PATCH_ID_INPUT and PATCH_ID_OUTPUT) before they were built from the
# camera configuration file:
# patch_sw_id of each position of the trigger traces, sectors in the order (1, 2, 3)
EXPECTED_PATCH_ID_INPUT = [
    204, 216, 180, 192, 229, 241, 205, 217, 254, 266, 230, 242,
    279, 291, 255, 267, 304, 316, 280, 292, 329, 341, 305, 317, 156, 168, 132,
    144, 181, 193, 157, 169, 206, 218, 182, 194, 231, 243, 207, 219, 256, 268,
    232, 244, 281, 293, 257, 269, 108, 120, 84, 96, 133, 145, 109, 121, 158,
    170, 134, 146, 183, 195, 159, 171, 208, 220, 184, 196, 233, 245, 209, 221,
    60, 72, 40, 50, 85, 97, 61, 73, 110, 122, 86, 98, 135, 147, 111, 123, 160,
    172, 136, 148, 185, 197, 161, 173, 24, 32, 12, 18, 41, 51, 25, 33, 62, 74,
    42, 52, 87, 99, 63, 75, 112, 124, 88, 100, 137, 149, 113, 125, 4, 8, 0, 2,
    13, 19, 5, 9, 26, 34, 14, 20, 43, 53, 27, 35, 64, 76, 44, 54, 89, 101, 65,
    77, 228, 239, 240, 252, 251, 262, 263, 275, 274, 285, 286, 298, 297, 308,
    309, 321, 320, 331, 332, 344, 343, 354, 355, 366, 253, 264, 265, 277, 276,
    287, 288, 300, 299, 310, 311, 323, 322, 333, 334, 346, 345, 356, 357, 368,
    367, 377, 378, 387, 278, 289, 290, 302, 301, 312, 313, 325, 324, 335, 336,
    348, 347, 358, 359, 370, 369, 379, 380, 389, 388, 396, 397, 404, 303, 314,
    315, 327, 326, 337, 338, 350, 349, 360, 361, 372, 371, 381, 382, 391, 390,
    398, 399, 406, 405, 411, 412, 417, 328, 339, 340, 352, 351, 362, 363, 374,
    373, 383, 384, 393, 392, 400, 401, 408, 407, 413, 414, 419, 418, 422, 423,
    426, 353, 364, 365, 376, 375, 385, 386, 395, 394, 402, 403, 410, 409, 415,
    416, 421, 420, 424, 425, 428, 427, 429, 430, 431, 215, 191, 227, 203, 167,
    143, 179, 155, 119, 95, 131, 107, 71, 49, 83, 59, 31, 17, 39, 23, 7, 1,
    11, 3, 238, 214, 250, 226, 190, 166, 202, 178, 142, 118, 154, 130, 94, 70,
    106, 82, 48, 30, 58, 38, 16, 6, 22, 10, 261, 237, 273, 249, 213, 189, 225,
    201, 165, 141, 177, 153, 117, 93, 129, 105, 69, 47, 81, 57, 29, 15, 37,
    21, 284, 260, 296, 272, 236, 212, 248, 224, 188, 164, 200, 176, 140, 116,
    152, 128, 92, 68, 104, 80, 46, 28, 56, 36, 307, 283, 319, 295, 259, 235,
    271, 247, 211, 187, 223, 199, 163, 139, 175, 151, 115, 91, 127, 103, 67,
    45, 79, 55, 330, 306, 342, 318, 282, 258, 294, 270, 234, 210, 246, 222,
    186, 162, 198, 174, 138, 114, 150, 126, 90, 66, 102, 78
]

EXPECTED_PATCH_ID_OUTPUT = [
    204, 216, 229, 241, 254, 266, 279, 291, 304, 316, 329,
    341, 180, 192, 205, 217, 230, 242, 255, 267, 280, 292, 305, 317, 156, 168,
    181, 193, 206, 218, 231, 243, 256, 268, 281, 293, 132, 144, 157, 169, 182,
    194, 207, 219, 232, 244, 257, 269, 108, 120, 133, 145, 158, 170, 183, 195,
    208, 220, 233, 245, 84, 96, 109, 121, 134, 146, 159, 171, 184, 196, 209,
    221, 60, 72, 85, 97, 110, 122, 135, 147, 160, 172, 185, 197, 40, 50, 61,
    73, 86, 98, 111, 123, 136, 148, 161, 173, 24, 32, 41, 51, 62, 74, 87, 99,
    112, 124, 137, 149, 12, 18, 25, 33, 42, 52, 63, 75, 88, 100, 113, 125, 4,
    8, 13, 19, 26, 34, 43, 53, 64, 76, 89, 101, 0, 2, 5, 9, 14, 20, 27, 35,
    44, 54, 65, 77, 228, 239, 251, 262, 274, 285, 297, 308, 320, 331, 343,
    354, 240, 252, 263, 275, 286, 298, 309, 321, 332, 344, 355, 366, 253, 264,
    276, 287, 299, 310, 322, 333, 345, 356, 367, 377, 265, 277, 288, 300, 311,
    323, 334, 346, 357, 368, 378, 387, 278, 289, 301, 312, 324, 335, 347, 358,
    369, 379, 388, 396, 290, 302, 313, 325, 336, 348, 359, 370, 380, 389, 397,
    404, 303, 314, 326, 337, 349, 360, 371, 381, 390, 398, 405, 411, 315, 327,
    338, 350, 361, 372, 382, 391, 399, 406, 412, 417, 328, 339, 351, 362, 373,
    383, 392, 400, 407, 413, 418, 422, 340, 352, 363, 374, 384, 393, 401, 408,
    414, 419, 423, 426, 353, 364, 375, 385, 394, 402, 409, 415, 420, 424, 427,
    429, 365, 376, 386, 395, 403, 410, 416, 421, 425, 428, 430, 431, 215, 191,
    167, 143, 119, 95, 71, 49, 31, 17, 7, 1, 227, 203, 179, 155, 131, 107, 83,
    59, 39, 23, 11, 3, 238, 214, 190, 166, 142, 118, 94, 70, 48, 30, 16, 6,
    250, 226, 202, 178, 154, 130, 106, 82, 58, 38, 22, 10, 261, 237, 213, 189,
    165, 141, 117, 93, 69, 47, 29, 15, 273, 249, 225, 201, 177, 153, 129, 105,
    81, 57, 37, 21, 284, 260, 236, 212, 188, 164, 140, 116, 92, 68, 46, 28,
    296, 272, 248, 224, 200, 176, 152, 128, 104, 80, 56, 36, 307, 283, 259,
    235, 211, 187, 163, 139, 115, 91, 67, 45, 319, 295, 271, 247, 223, 199,
    175, 151, 127, 103, 79, 55, 330, 306, 282, 258, 234, 210, 186, 162, 138,
    114, 90, 66, 342, 318, 294, 270, 246, 222, 198, 174, 150, 126, 102, 78
]


def test_patch_ids_from_camera_config():
    patches = read_trigger_patches()
    assert set(patches) == set(CAMERA_CONFIG_TRIGGER_COLUMNS)
    assert all(len(values) == N_PATCHES for values in patches.values())
    assert np.array_equal(patches["patch_sw_id"], np.arange(N_PATCHES))
    # 3 sectors of 36 modules of 4 patches
    assert sorted(np.unique(patches["micro_crate"], return_counts=True)[1]) == [PATCHES_PER_SECTOR] * 3
    assert set(patches["patch_in_mod_fw"]) == {0, 1, 2, 3}
    assert set(patches["patch_trig_out"]) == set(range(PATCHES_PER_SECTOR))

    patch_ids = build_patch_ids(patches)
    assert patch_ids["input"] == EXPECTED_PATCH_ID_INPUT
    assert patch_ids["output"] == EXPECTED_PATCH_ID_OUTPUT
    assert build_patch_ids() == patch_ids


def test_build_patch_ids_refuses_two_patches_at_the_same_position():
    patches = read_trigger_patches()
    patches["patch_trig_out"][1] = patches["patch_trig_out"][0]
    patches["micro_crate"][1] = patches["micro_crate"][0]
    with pytest.raises(ValueError, match="same position"):
        build_patch_ids(patches)


def test_read_trigger_patches_refuses_inconsistent_patches(tmp_path):
    with open(CAMERA_CONFIG_FILE) as f:
        lines = f.readlines()
    first_pixel = next(i for i, line in enumerate(lines) if line.strip() and not line.startswith("#"))
    row = lines[first_pixel].split()
    row[CAMERA_CONFIG_TRIGGER_COLUMNS["patch_trig_out"]] = "143"  # not the value of the 2 other pixels
    lines[first_pixel] = " ".join(row) + "\n"
    config_file = tmp_path / "camera_config.cfg"
    config_file.write_text("".join(lines))
    with pytest.raises(ValueError, match="different readout values"):
        read_trigger_patches(config_file)


def test_geometry_in_engineering_frame():
    """the CTS geometry is in the EngineeringCameraFrame: in the CameraFrame, it is the one of the subarray"""
    import astropy.units as u
    from ctapipe.coordinates import CameraFrame, EngineeringCameraFrame
    from ctapipe.instrument import SubarrayDescription

    from sst1mpipe.resources import SUBARRAY_FILE

    assert isinstance(GEOMETRY.frame, EngineeringCameraFrame)
    in_camera_frame = GEOMETRY.transform_to(CameraFrame())
    for tel in SubarrayDescription.from_hdf(SUBARRAY_FILE).tel.values():
        subarray_geometry = tel.camera.geometry
        assert isinstance(subarray_geometry.frame, CameraFrame)
        np.testing.assert_allclose(in_camera_frame.pix_x.to_value(u.mm), subarray_geometry.pix_x.to_value(u.mm), atol=1e-3)
        np.testing.assert_allclose(in_camera_frame.pix_y.to_value(u.mm), subarray_geometry.pix_y.to_value(u.mm), atol=1e-3)
