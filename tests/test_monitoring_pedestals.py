import numpy as np
import pytest
from ctapipe.containers import MonitoringCameraContainer, PedestalContainer
from ctapipe.image import ImageProcessor

from sst1mpipe.io import load_config
from sst1mpipe.io.containers import SST1MArrayEventContainer
from sst1mpipe.utils import get_subarray
from sst1mpipe.utils.monitoring_pedestals import (
    DL1PedestalMonitor,
    R0PedestalMonitor,
    SlidingWindowMonitor,
)

SUBARRAY = get_subarray()
TEL_ID = 21
N_PIXELS = SUBARRAY.tel[TEL_ID].camera.geometry.n_pixels
N_SAMPLES = 50


def make_event(rng, time_s, tel_id=TEL_ID, image_std=1.0):
    # NOTE: all SST1MArrayEventContainer share the same sst1m container
    event = SST1MArrayEventContainer()
    r0 = event.r0.tel[tel_id]
    r0.waveform = rng.normal(300, 5, (N_PIXELS, N_SAMPLES))
    r0.event_time = int(time_s * 1e9)
    event.dl1.tel[tel_id].image = rng.normal(0, image_std, N_PIXELS)
    return event


class ImageMonitor(SlidingWindowMonitor):
    """Minimal subclass monitoring the images in event.mon.tel[tel_id].pedestal"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.n_computations = 0

    def add_event(self, event, tel_id):
        self._append(event, tel_id, event.dl1.tel[tel_id].image.copy())

    def _container(self, event, tel_id):
        return event.mon.tel[tel_id].pedestal

    def _compute_statistics(self, values):
        self.n_computations += 1
        return values.mean(axis=0), np.median(values, axis=0), values.std(axis=0)


# ---------------------------------------------------------------------------
# SlidingWindowMonitor
# ---------------------------------------------------------------------------


def test_base_class_is_abstract():
    with pytest.raises(TypeError):
        SlidingWindowMonitor(subarray=SUBARRAY)

    class MissingContainer(SlidingWindowMonitor):
        def add_event(self, event, tel_id):
            pass

        def _compute_statistics(self, values):
            pass

    with pytest.raises(TypeError):
        MissingContainer(subarray=SUBARRAY)


def test_sliding_window_per_telescope():
    rng = np.random.default_rng(0)
    monitor = ImageMonitor(subarray=SUBARRAY, n_events=[("id", 21, 3), ("id", 22, 5)])

    for i in range(8):
        monitor.add_event(make_event(rng, i, tel_id=21), 21)
    for i in range(2):
        monitor.add_event(make_event(rng, i, tel_id=22), 22)

    assert monitor.n_buffered(21) == 3 and monitor.processed_events[21] == 8
    assert monitor.n_buffered(22) == 2 and monitor.processed_events[22] == 2
    assert monitor.n_buffered(1) == 0


def test_call_adds_event_and_fills_container():
    rng = np.random.default_rng(1)
    monitor = ImageMonitor(subarray=SUBARRAY, n_events=3)
    images = []
    for i in range(5):
        event = make_event(rng, time_s=10 + i)
        images.append(event.dl1.tel[TEL_ID].image)
        monitor(event, TEL_ID)

    container = event.mon.tel[TEL_ID].pedestal
    images = np.array(images[-3:])
    assert container.n_events == 3
    assert container.sample_time.to_value("s") == 13
    assert container.sample_time_min.to_value("s") == 12
    assert container.sample_time_max.to_value("s") == 14
    np.testing.assert_allclose(container.charge_mean, images.mean(axis=0))
    np.testing.assert_allclose(container.charge_median, np.median(images, axis=0))
    np.testing.assert_allclose(container.charge_std, images.std(axis=0))


def test_empty_window_does_not_fill_container():
    event = SST1MArrayEventContainer()
    ImageMonitor(subarray=SUBARRAY).fill_monitoring(event, TEL_ID)

    assert event.mon.tel[TEL_ID].pedestal.n_events == -1
    assert event.mon.tel[TEL_ID].pedestal.charge_std is None


def test_statistics_only_computed_for_new_events():
    rng = np.random.default_rng(2)
    monitor = ImageMonitor(subarray=SUBARRAY)
    event = make_event(rng, time_s=0)

    monitor(event, TEL_ID)
    monitor.fill_monitoring(event, TEL_ID)
    monitor.fill_monitoring(event, TEL_ID)
    assert monitor.n_computations == 1

    monitor(make_event(rng, time_s=1), TEL_ID)
    assert monitor.n_computations == 2


# ---------------------------------------------------------------------------
# R0PedestalMonitor and DL1PedestalMonitor
# ---------------------------------------------------------------------------


def test_monitoring_containers_are_ctapipe_compatible():
    mon = SST1MArrayEventContainer().mon.tel[TEL_ID]
    assert isinstance(mon, MonitoringCameraContainer)
    assert isinstance(mon.r0, PedestalContainer)


def test_configuration_from_data_config():
    config = load_config(None, ismc=False)
    assert R0PedestalMonitor(subarray=SUBARRAY, config=config).n_events.tel[TEL_ID] == 100
    assert DL1PedestalMonitor(subarray=SUBARRAY, config=config).n_events.tel[TEL_ID] == 1000


def test_r0_pedestal_monitor():
    rng = np.random.default_rng(3)
    monitor = R0PedestalMonitor(subarray=SUBARRAY, n_events=3)
    samples = []
    for i in range(5):
        event = make_event(rng, time_s=i)
        samples.append(event.r0.tel[TEL_ID].waveform)
        monitor(event, TEL_ID)

    container = event.mon.tel[TEL_ID].r0
    means = np.array(samples[-3:]).mean(axis=2)
    stds = np.array(samples[-3:]).std(axis=2)
    np.testing.assert_allclose(container.charge_mean, means.mean(axis=0))
    np.testing.assert_allclose(container.charge_median, np.median(means, axis=0))
    np.testing.assert_allclose(container.charge_std, stds.mean(axis=0))
    # the dl1 container is not filled
    assert event.mon.tel[TEL_ID].pedestal.charge_std is None


def test_r0_pedestal_monitor_ignores_masked_pixels():
    rng = np.random.default_rng(4)
    monitor = R0PedestalMonitor(subarray=SUBARRAY)
    mask = np.zeros(N_PIXELS, dtype=bool)
    mask[:10] = True

    event = make_event(rng, time_s=0)
    samples = event.r0.tel[TEL_ID].waveform.copy()
    monitor(event, TEL_ID, cleaning_mask=mask)
    # the event is not modified
    np.testing.assert_array_equal(event.r0.tel[TEL_ID].waveform, samples)

    event = make_event(rng, time_s=1)
    second = event.r0.tel[TEL_ID].waveform[mask]
    monitor(event, TEL_ID)

    # the masked pixels only use the second event
    container = event.mon.tel[TEL_ID].r0
    np.testing.assert_allclose(container.charge_mean[mask], second.mean(axis=1))
    np.testing.assert_allclose(container.charge_std[mask], second.std(axis=1))


def test_dl1_pedestal_monitor_raises_cleaning_threshold():
    """NSBImageCleaner uses the std of the pedestal images in event.mon.tel[tel_id].pedestal"""
    rng = np.random.default_rng(5)
    image_processor = ImageProcessor(subarray=SUBARRAY, config=load_config(None, ismc=False))
    monitor = DL1PedestalMonitor(subarray=SUBARRAY)

    image = np.full(N_PIXELS, 12.0)
    times = np.zeros(N_PIXELS)
    event = SST1MArrayEventContainer()
    assert image_processor.clean(TEL_ID, image, arrival_times=times, monitoring=event.mon.tel[TEL_ID]).all()

    for i in range(50):
        event = make_event(rng, time_s=i, image_std=10)
        monitor(event, TEL_ID)

    # 2.5 * std ~ 25 p.e. > 12 p.e.
    assert event.mon.tel[TEL_ID].r0.charge_std is None
    assert not image_processor.clean(TEL_ID, image, arrival_times=times, monitoring=event.mon.tel[TEL_ID]).any()
