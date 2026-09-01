from importlib.resources import files
import astropy.units as u

import numpy as np
from ctapipe.image import ImageProcessor, dilate
from ctapipe.instrument import SubarrayDescription
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist

from sklearn.cluster import DBSCAN

from ctapipe.image.cleaning import apply_time_delta_cleaning, tailcuts_clean, ImageCleaner
from ctapipe.image.toymodel import Gaussian

from traitlets.config import Config

from sst1mpipe.instrument.camera import Camera
from sst1mpipe.utils.cleaning import DBSCANImageCleaner, TimeDBSCANImageCleaner, DBSCANImageCleaner3D, clean_dbscan_fast

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

def test_fast_dbscan():

    max_pe = 100
    epsilon = 40

    x = np.column_stack([GEOMETRY.pix_x.value, GEOMETRY.pix_y.value])
    d = cdist(x, x)
    connectivity = d <= epsilon

    for _ in range(100):

        image = np.random.randint(0, max_pe, GEOMETRY.n_pixels)
        dbscan = DBSCAN(eps=epsilon, min_samples=max_pe // 2 * 10, metric='precomputed').fit(d, sample_weight=image)
        mask_dbscan = (dbscan.labels_ >= 0)
        mask = clean_dbscan_fast(connectivity, weights=image, min_points=max_pe // 2 * 10)

        assert (mask == mask_dbscan).all()

def test_dbscan_similar_to_tailcuts():
    """This test checks that the dbscan cleaning is identical to the
    tailcuts cleaning under the condition that epsilon = 0"""

    image, signal, noise = Gaussian(x=0 * u.mm, y=0 * u.mm, length=40 * u.mm, width=10 * u.mm, psi=0 * u.rad).generate_image(camera=CAMERA.geometry)
    x = np.column_stack([GEOMETRY.pix_x.value, GEOMETRY.pix_y.value])
    d = cdist(x, x)
    connectivity = d == 0.0

    for threshold in range(10):

        mask_tailcuts = tailcuts_clean(CAMERA.geometry, image, picture_thresh=threshold,
                                       boundary_thresh=-np.inf, keep_isolated_pixels=True,
                                       min_number_picture_neighbors=0)
        mask_dbscan = clean_dbscan_fast(connectivity, weights=image, min_points=threshold)
        mask_dbscan = dilate(CAMERA.geometry, mask_dbscan)

        assert (mask_dbscan == mask_tailcuts).all()
