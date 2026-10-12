import numpy as np
from sst1mpipe.instrument.camera import Camera

CAMERA = Camera()
GEOMETRY = CAMERA.geometry

def test_geometry_pixel_has_units():

    assert GEOMETRY.pix_x.unit.physical_type == "length"
    assert GEOMETRY.pix_y.unit.physical_type == "length"
    assert GEOMETRY.pix_area.unit.physical_type == "area"

def test_get_neighbors_sparse_matrix():

    assert isinstance(GEOMETRY.neighbor_matrix, np.ndarray)
