
import numpy as np
import pytest
from ctapipe.containers import MonitoringCameraContainer
from ctapipe.image import ImageProcessor, apply_time_delta_cleaning, tailcuts_clean
from ctapipe.image.cleaning import NSBImageCleaner
from ctapipe.instrument import SubarrayDescription

from sst1mpipe.io import load_config
from sst1mpipe.resources import SUBARRAY_FILE

SUBARRAY = SubarrayDescription.from_hdf(SUBARRAY_FILE)
TEL_ID = 21
GEOMETRY = SUBARRAY.tel[TEL_ID].camera.geometry


def make_cleaner():
    config = load_config(None, ismc=False)
    return ImageProcessor(subarray=SUBARRAY, config=config).clean


def two_islands_event():
    """Two bright islands in a quiet camera, all pixels in time."""
    image = np.zeros(GEOMETRY.n_pixels)
    big = np.array([0, *GEOMETRY.neighbors[0]])
    far = np.argmax(np.hypot(GEOMETRY.pix_x - GEOMETRY.pix_x[0], GEOMETRY.pix_y - GEOMETRY.pix_y[0]))
    small = np.array([far, *GEOMETRY.neighbors[far][:2]])
    image[big] = 50
    image[small] = 50
    return image, np.full(GEOMETRY.n_pixels, 20.0), big, small


@pytest.mark.parametrize("ismc", [False, True])
def test_default_configs_use_nsb_image_cleaner(ismc):
    config = load_config(None, ismc=ismc)
    cleaner = ImageProcessor(subarray=SUBARRAY, config=config).clean

    assert isinstance(cleaner, NSBImageCleaner)
    assert cleaner.bright_cleaning_threshold.tel[TEL_ID] is None


def test_without_pedestal_std_is_tailcuts_and_time_delta():
    cleaner = make_cleaner()
    rng = np.random.default_rng(0)
    image = rng.exponential(3, GEOMETRY.n_pixels)
    image[:200] += 20
    times = rng.normal(20, 5, GEOMETRY.n_pixels)

    expected = tailcuts_clean(GEOMETRY, image, picture_thresh=8, boundary_thresh=4, min_number_picture_neighbors=2)
    expected = apply_time_delta_cleaning(GEOMETRY, expected, times, min_number_neighbors=1, time_limit=8)

    for monitoring in (None, MonitoringCameraContainer()):
        mask = cleaner(TEL_ID, image, arrival_times=times, monitoring=monitoring)
        assert expected.any()
        np.testing.assert_array_equal(mask, expected)


def test_pedestal_std_raises_picture_threshold():
    cleaner = make_cleaner()
    image, times, big, small = two_islands_event()

    monitoring = MonitoringCameraContainer()
    pedestal_std = np.zeros(GEOMETRY.n_pixels)
    pedestal_std[small] = 25  # 2.5 * 25 > 50 p.e.
    monitoring.pedestal.charge_std = pedestal_std

    mask = cleaner(TEL_ID, image, arrival_times=times)
    assert mask[big].all() and mask[small].all()

    mask = cleaner(TEL_ID, image, arrival_times=times, monitoring=monitoring)
    assert mask[big].all()
    assert not mask[small].any()
