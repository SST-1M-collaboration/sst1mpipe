import logging
from abc import abstractmethod
from collections import deque
from copy import deepcopy

import astropy.units as u
from astropy.time import Time
import numpy as np
from ctapipe.calib import CameraCalibrator
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import IntTelescopeParameter
from ctapipe.image import ImageProcessor

from sst1mpipe.calib import (
    R0R1Calibrator,
    ImageSaturationCorrector,
)
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils import get_subarray

MON_EVT_TYPE = 8
MASKED_VALUE = -100


class SlidingWindowMonitor(TelescopeComponent):
    """
    Base class keeping per pixel quantities of the last ``n_events`` events of each
    telescope in a sliding window, and filling a monitoring container of
    ``event.mon.tel[tel_id]`` with their statistics.

    Subclasses define which quantities are extracted from an event (`add_event`),
    how their statistics are computed (`_compute_statistics`) and which container
    is filled (`_container`).
    """

    n_events = IntTelescopeParameter(
        default_value=100, help="Number of events in the sliding window"
    ).tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        self._timestamps = {}
        self._values = {}
        self._statistics = {}
        self.processed_events = {}

    def __call__(self, event, tel_id, **kwargs):
        """
        Add the event to the sliding window and fill the monitoring container.
        ``kwargs`` are passed to `add_event`.
        """
        self.add_event(event, tel_id, **kwargs)
        self.fill_monitoring(event, tel_id)

    @abstractmethod
    def add_event(self, event, tel_id, **kwargs):
        """
        Add an event to the sliding window, extracting its quantities and
        passing them to `_append`.
        """

    @abstractmethod
    def _container(self, event, tel_id):
        """The container of ``event`` to fill"""

    @abstractmethod
    def _compute_statistics(self, values):
        """
        Return charge_mean, charge_median, charge_std of ``values``, an array
        of shape (n_buffered, ...) of the quantities given to `_append`
        """

    def n_buffered(self, tel_id):
        """Number of events in the sliding window"""
        return len(self._timestamps.get(tel_id, ()))

    def _append(self, event, tel_id, values):
        if tel_id not in self._timestamps:
            self._timestamps[tel_id] = deque(maxlen=self.n_events.tel[tel_id])
            self._values[tel_id] = deque(maxlen=self.n_events.tel[tel_id])
            self.processed_events[tel_id] = 0

        self._timestamps[tel_id].append(event.r0.tel[tel_id].event_time)
        self._values[tel_id].append(values)
        self.processed_events[tel_id] += 1
        self._statistics.pop(tel_id, None)

    def fill_monitoring(self, event, tel_id):
        """
        Fill the monitoring container of ``event.mon.tel[tel_id]``.
        Nothing is done if the sliding window is empty.
        """
        if self.n_buffered(tel_id) == 0:
            return

        # statistics are only recomputed when a new event is added
        if tel_id not in self._statistics:
            self._statistics[tel_id] = self._compute_statistics(np.array(self._values[tel_id]))
        charge_mean, charge_median, charge_std = self._statistics[tel_id]

        timestamps = self._timestamps[tel_id]
        container = self._container(event, tel_id)
        container.n_events = len(timestamps)
        # the ctapipe PedestalContainer stores the times as Quantity [s] (unix TAI)
        container.sample_time = Time(timestamps).mean().unix_tai * u.s
        container.sample_time_min = timestamps[0].unix_tai * u.s
        container.sample_time_max = timestamps[-1].unix_tai * u.s
        container.charge_mean = charge_mean
        container.charge_median = charge_median
        container.charge_std = charge_std


class R0PedestalMonitor(SlidingWindowMonitor):
    """
    Statistics of the ADC samples of the pedestal events, filled in
    ``event.mon.tel[tel_id].r0``. Used for the voltage drop
    correction and the identification of dead pixels.
    """

    def add_event(self, event, tel_id, cleaning_mask=None):
        """
        Add the ADC samples of a pedestal event. Pixels in ``cleaning_mask``
        (e.g. Cherenkov pixels of a fake pedestal) are not used.
        """
        samples = event.r0.tel[tel_id].waveform
        if cleaning_mask is not None:
            samples = samples.astype(np.float64)
            samples[cleaning_mask] = MASKED_VALUE
        self._append(event, tel_id, np.stack([samples.mean(axis=-1), samples.std(axis=-1)]))

    def _container(self, event, tel_id):
        return event.mon.tel[tel_id].r0

    def _compute_statistics(self, values):
        means = np.ma.masked_values(values[:, 0], MASKED_VALUE)
        # masked pixels have a null standard deviation
        stds = np.ma.masked_values(values[:, 1], 0)
        return (
            means.mean(axis=0).filled(np.nan),
            np.ma.median(means, axis=0).filled(np.nan),
            stds.mean(axis=0).filled(0),
        )


class DL1PedestalMonitor(SlidingWindowMonitor):
    """
    Statistics of the calibrated images (p.e.) of the pedestal events, filled in
    ``event.mon.tel[tel_id].pedestal``. Used by `ctapipe.image.cleaning.NSBImageCleaner`.
    """

    n_events = IntTelescopeParameter(
        default_value=1000, help="Number of events in the sliding window"
    ).tag(config=True)

    def add_event(self, event, tel_id):
        """Add the calibrated image of a pedestal event."""
        self._append(event, tel_id, np.array(event.dl1.tel[tel_id].image, dtype=np.float64))

    def _container(self, event, tel_id):
        return event.mon.tel[tel_id].pedestal

    def _compute_statistics(self, values):
        return values.mean(axis=0), np.median(values, axis=0), values.std(axis=0)


def load_first_pedestals(r0_monitor, dl1_monitor, input_file, config, max_events=100000):
    """
    Fill the monitors with the first pedestal events of ``input_file``.
    The ADC samples of the first pedestal events are needed to estimate the voltage
    drop, which is used to calibrate the images of the pedestal events.

    If there are no pedestal events, shower/NSB events with their Cherenkov pixels
    masked out (fake pedestals) are used instead.

    Returns
    -------
    pedestals_in_file: bool
        False if fake pedestals were used
    """
    source = SST1MEventSource(input_url=input_file, config=config, max_events=max_events)
    source._subarray = get_subarray()
    tel = None

    for event in source:
        tel = event.trigger.tels_with_trigger[0]
        if event.r0.tel[tel]._event_type.value == MON_EVT_TYPE:
            r0_monitor.add_event(event, tel)
        if r0_monitor.n_buffered(tel) >= r0_monitor.n_events.tel[tel]:
            break

    if tel is not None and r0_monitor.n_buffered(tel) > 0:
        _load_first_images(r0_monitor, dl1_monitor, source, tel, config)
        logging.info("%d pedestal events loaded in buffer", r0_monitor.n_buffered(tel))
        return True

    logging.warning("No pedestal events found in firsts events. Cleaned shower/NSB events used instead.")
    tel = _load_first_fake_pedestals(r0_monitor, dl1_monitor, input_file, config)
    if tel is not None:
        logging.info("%d fake pedestal events loaded in buffer", r0_monitor.n_buffered(tel))
    return False


def _load_first_images(r0_monitor, dl1_monitor, source, tel, config):

    r1_dl1_calibrator = CameraCalibrator(subarray=source.subarray, config=config)
    calibrator_r0_r1 = R0R1Calibrator(subarray=source.subarray, config=config)
    image_saturation_corrector = ImageSaturationCorrector(subarray=source.subarray, config=config)

    for event in source:

        if event.r0.tel[tel]._event_type.value != MON_EVT_TYPE:
            continue

        # here we apply gain drop correction
        r0_monitor.fill_monitoring(event, tel)
        calibrator_r0_r1(event, tel)
        r1_dl1_calibrator(event)

        # Integration correction of saturated pixels
        image_saturation_corrector(event, tel)

        dl1_monitor.add_event(event, tel)
        if dl1_monitor.n_buffered(tel) >= dl1_monitor.n_events.tel[tel]:
            break


def _load_first_fake_pedestals(r0_monitor, dl1_monitor, input_file, config, max_events=10):

    # Here (for the first few events) we use just the simple ImageProcessor, nothing fancy
    config = deepcopy(config)
    config["ImageProcessor"]["image_cleaner_type"] = "TailcutsImageCleaner"

    source = SST1MEventSource(input_url=input_file, config=config, max_events=max_events)
    source._subarray = get_subarray()
    r1_dl1_calibrator = CameraCalibrator(subarray=source.subarray, config=config)
    image_processor = ImageProcessor(subarray=source.subarray, config=config)
    calibrator_r0_r1 = R0R1Calibrator(subarray=source.subarray, config=config)
    image_saturation_corrector = ImageSaturationCorrector(subarray=source.subarray, config=config)
    tel = None

    def clean(event):
        calibrator_r0_r1(event, tel)
        r1_dl1_calibrator(event)
        image_processor(event)
        return event.dl1.tel[tel].image_mask

    for event in source:
        if tel is None:
            tel = event.trigger.tels_with_trigger[0]

        cleaning_mask = clean(event)
        # Arbitrary cut, just to prevent too big showers from being used
        if sum(cleaning_mask) < 20:
            r0_monitor.add_event(event, tel, cleaning_mask=cleaning_mask)

    if tel is None or r0_monitor.n_buffered(tel) == 0:
        return tel

    # to treat images we need to estimate gain drop from pedestals, so calibrate the events once more
    for event in source:

        r0_monitor.fill_monitoring(event, tel)
        cleaning_mask = clean(event)
        if sum(cleaning_mask) < 20:
            # Integration correction of saturated pixels - done only here because the fake pedestals must match in both loops
            image_saturation_corrector(event, tel)
            dl1_monitor.add_event(event, tel)
            if dl1_monitor.n_buffered(tel) >= dl1_monitor.n_events.tel[tel]:
                break

    return tel
