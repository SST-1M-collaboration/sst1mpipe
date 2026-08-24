import numpy as np
from ctapipe.image import ImageProcessor
from scipy.sparse import csr_matrix

from ctapipe.image.cleaning import apply_time_delta_cleaning, ImageCleaner
from traitlets.config import Config

from sst1mpipe.instrument.camera import Camera
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils.cleaning import DBSCANImageCleaner, TimeDBSCANImageCleaner, DBSCANImageCleaner3D

CAMERA = Camera()
GEOMETRY = CAMERA.geometry
SUBARRAY = SST1MEventSource.create_subarray()

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
                ["id", 1, 38],
                ["id", 2, 38]
            ],
            "minimum_pe": [
                ["id", 1, 30],
                ["id", 2, 30]
            ]
        }}})

    ImageProcessor(subarray=SUBARRAY, config=config)
