from importlib.resources import files
import astropy.units as u

import numpy as np
from ctapipe.image import ImageProcessor
from ctapipe.instrument import SubarrayDescription
from scipy.sparse import csr_matrix

from ctapipe.image.cleaning import apply_time_delta_cleaning, ImageCleaner
from traitlets.config import Config

from sst1mpipe.instrument.camera import Camera
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils.cleaning import DBSCANImageCleaner, TimeDBSCANImageCleaner, DBSCANImageCleaner3D

SUBARRAY_FILE = files('sst1mpipe.data').joinpath(
    'sst1m_array.h5'
)

CAMERA = Camera()
GEOMETRY = CAMERA.geometry
SUBARRAY =  SubarrayDescription.from_hdf(SUBARRAY_FILE)

def test_sparse_matrix():

    a = csr_matrix([1])

    assert a.A == a.toarray()

def test_apply_time_delta_cleaning():

    mask = np.ones(GEOMETRY.n_pixels, dtype=bool)
    times = np.ones(GEOMETRY.n_pixels)

    apply_time_delta_cleaning(GEOMETRY, mask, times, min_number_neighbors = 1, time_limit = 8)


def test_load_DBSCANCleaning_from_name():


    cleaner = ImageCleaner.from_name(subarray=SUBARRAY, name='DBSCANImageCleaner')
    assert isinstance(cleaner, DBSCANImageCleaner)
    cleaner = ImageCleaner.from_name(subarray=SUBARRAY, name='TimeDBSCANImageCleaner')
    assert isinstance(cleaner, TimeDBSCANImageCleaner)
    cleaner = ImageCleaner.from_name(subarray=SUBARRAY, name='DBSCANImageCleaner3D')
    assert isinstance(cleaner, DBSCANImageCleaner3D)

def test_DBSCAN_configurable_from_image_processor():

    config = Config({"ImageProcessor": {
        "image_cleaner_type": "DBSCANImageCleaner",
        "use_telescope_frame": False,
        "DBSCANImageCleaner": {
            "epsilon_r": [
                ["id", 21, 38],
                ["id", 22, 38]
            ],
            "minimum_pe": [
                ["id", 21, 30],
                ["id", 22, 30]
            ]
        }}})

    ImageProcessor(subarray=SUBARRAY, config=config)

def test_DBSCANImageCleaner():

    image = np.ones(GEOMETRY.n_pixels)

    epsilon_r = (np.sqrt(GEOMETRY.pix_area[0] / np.sqrt(3.0) / 2)).to(u.mm).value * 2 # 2 pixel to pixel distance


    n_pixels = [GEOMETRY.n_pixels, 1284, 1282, 1242, 857, 345, 133, 0, 0]

    for i in range(9):

        cleaner = DBSCANImageCleaner(subarray=SUBARRAY, minimum_pe=[("id", 21, i+1), ("id", 22, i+1)],
                                     epsilon_r=[("id", 21, epsilon_r), ("id", 22, epsilon_r)])
        mask = cleaner(21, image)
        assert mask.sum() == n_pixels[i]
