import numpy as np
from scipy.sparse import csr_matrix

from ctapipe.image.cleaning import apply_time_delta_cleaning

from sst1mpipe.instrument.camera import Camera

CAMERA = Camera()
GEOMETRY = CAMERA.geometry

def test_sparse_matrix():

    a = csr_matrix([1])

    assert a.A == a.toarray()

def test_apply_time_delta_cleaning():

    mask = np.ones(GEOMETRY.n_pixels, dtype=bool)
    times = np.ones(GEOMETRY.n_pixels)

    apply_time_delta_cleaning(GEOMETRY, mask, times, min_number_neighbors = 1, time_limit = 8)

