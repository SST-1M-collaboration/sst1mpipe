import logging

import astropy.units as u
import numpy as np
import pandas as pd

from ctapipe.containers import PixelStatus, R1CameraContainer
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import (
    BoolTelescopeParameter,
    CaselessStrEnum,
    Float,
    FloatTelescopeParameter,
    Int,
    IntTelescopeParameter,
    List,
    Path,
    TelescopeParameter,
)

from sst1mpipe.utils import VAR_to_Idrop
from sst1mpipe.resources import CALIBRATION_DIR, WINDOW_DIR



# Transmittance correction factors of the camera windows, measured in the lab
DEFAULT_WINDOW_FILES = {
    21: 'corr_factor_1st_wdw.txt',
    22: 'corr_factor_2nd_wdw.txt',
}


def read_window_transmittance(window_file):
    """Correction factor of the window transmittance of each pixel (file with pixel_id, correction_factor)"""
    return np.loadtxt(window_file, unpack=True, skiprows=1, usecols=1)


def get_default_window(telescope=None):
    """
    Provides default window transmissivity file,
    used in the case when it is not defined in the
    configuration file.

    Parameters
    ----------
    telescope: int
        Telescope number as in
        event.trigger.tels_with_trigger (1, 2 for the simulations)

    Returns
    -------
    window_corr: numpy.ndarray
    window_file: string

    """
    window_file = (WINDOW_DIR / DEFAULT_WINDOW_FILES[{1: 21, 2: 22}.get(telescope, telescope)])
    logging.info('Window file used: %s', window_file)
    return read_window_transmittance(window_file), window_file


class ImageSaturationCorrector(TelescopeComponent):
    r"""
    Correction of the charges and peak times of the saturated pixels, applied after the
    image extraction (`~ctapipe.calib.CameraCalibrator`).

    The saturated pixels are the pixels flagged with `~ctapipe.containers.PixelStatus.SATURATED`
    in the pixel status of R1 (by the `R0R1Calibrator`). The standard integration window does
    not perform well for these pulses: their charge is the sum of the R1 waveform (p.e.) from the
    first sample above ``integration_level`` times the maximum, to the first sample below it after
    the maximum (or the end of the readout window). Their peak time is the middle of the ADC
    samples (``event.r0``, pedestal subtracted) above ``peak_time_level``.
    """

    integration_level = FloatTelescopeParameter(
        default_value=0.2,
        help="Fraction of the maximum defining the integration window of the saturated pulses",
    ).tag(config=True)

    peak_time_level = FloatTelescopeParameter(
        default_value=2500.0,
        help=(
            "ADC (pedestal subtracted): the peak time of a saturated pulse is the middle"
            " of its R0 samples above this level"
        ),
    ).tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        self.n_saturated_events = {}

    def __call__(self, event, tel_id=None):
        """
        Correct the saturated pixels of the telescope ``tel_id``, or of all the telescopes
        with DL1 data. Returns True if saturated pixels were corrected.
        """
        tel_ids = list(event.dl1.tel.keys()) if tel_id is None else [tel_id]
        saturated = False
        for tel in tel_ids:
            if self._correct_telescope(event, tel):
                self.n_saturated_events[tel] = self.n_saturated_events.get(tel, 0) + 1
                saturated = True
        return saturated

    @staticmethod
    def saturated_pixels(event, tel_id):
        """Pixels flagged as saturated in the pixel status of R1"""
        pixel_status = event.r1.tel[tel_id].pixel_status
        if pixel_status is None:
            return np.array([], dtype=int)
        return np.flatnonzero(pixel_status & PixelStatus.SATURATED)

    def _correct_telescope(self, event, tel_id):
        saturated = self.saturated_pixels(event, tel_id)
        if len(saturated) == 0:
            return False

        integration_level = self.integration_level.tel[tel_id]
        peak_time_level = self.peak_time_level.tel[tel_id]
        sample_time = (1 / self.subarray.tel[tel_id].camera.readout.sampling_rate).to_value(u.ns)

        dl1 = event.dl1.tel[tel_id]
        r0 = event.r0.tel[tel_id]
        r1_waveform = event.r1.tel[tel_id].waveform[0]
        for pixel in saturated:
            # the R1 waveform is proportional to the ADC samples of the pixel: same window
            w = r1_waveform[pixel]
            max_charge = w.max()
            index_max = np.argmax(w)
            integration_start = np.flatnonzero(w >= integration_level * max_charge)[0]
            # first drop below the level after the maximum, to avoid secondary peaks.
            # If there is none, the pulse extends over the readout window
            after_max = np.flatnonzero((w < integration_level * max_charge) & (np.arange(len(w)) > index_max))
            integration_stop = after_max[0] if len(after_max) > 0 else len(w) - 1
            dl1.image[pixel] = w[integration_start:integration_stop + 1].sum()

            adc = r0.waveform[0, pixel] - r0.pedestal[pixel]
            above_level = np.flatnonzero(adc > peak_time_level)
            dl1.peak_time[pixel] = (len(above_level) / 2 + above_level[0]) * sample_time
        return True


# Calibration parameters averaged from all darks taken between
# March 2023 and June 2024 (TEL1) and Sep 2023 and June 2024 (TEL2).
# Based on TT's study, there is a relative difference in the dc_to_pe
# factor between individual darks on the level of 5%, showing a
# slowly decreasing trend. In the future, we may start producing
# calibration files "per season" to mitigate the systematic uncertainty,
# but per-night is not necessary. TT also confirms that dc_to_pe
# does not depend on the level of DCR (the camera temperature).
DEFAULT_CALIBRATION_FILES = {
    21: 'averaged_calib_param_v2_2023_2024_tel1.h5',
    22: 'averaged_calib_param_v2_2023_2024_tel2.h5',
}
DEFAULT_CALIBRATION_FILES[1] = DEFAULT_CALIBRATION_FILES[21]
DEFAULT_CALIBRATION_FILES[2] = DEFAULT_CALIBRATION_FILES[22]

VOLTAGE_DROP_CORRECTIONS = ("none", "global", "pixelwise")


class R0R1Calibrator(TelescopeComponent):
    """
    R0 -> R1 calibration of the SST-1M telescopes. For the observed data, fills ``event.r1.tel[tel_id]``
    (`~ctapipe.containers.R1CameraContainer`, waveform in p.e. of shape
    (n_channels=1, n_pixels, n_samples)) from ``event.r0.tel[tel_id]``:

    1. subtraction of the pedestal computed by DigiCam (``r0.pedestal``)
    2. conversion from ADC to p.e. with the ``dc_to_pe`` of the calibration file
    3. voltage drop correction (``voltage_drop_correction``), from the std of the ADC
       samples of the pedestal events in ``event.mon.tel[tel_id].r0``
       (see `sst1mpipe.utils.monitoring_pedestals.R0PedestalMonitor`)
    4. window transmittance correction: division by the correction factor of each pixel
       of ``window_transmittance_file``
    5. bad pixels: pixels with bad calibration parameters (``flag_bad_calibration_pixels``),
       dead pixels (``flag_dead_pixels``) and the ``bad_pixels`` are set to 0 in the R1 waveforms and
       flagged in ``event.mon.tel[tel_id].pixel_status``, so that their charge is
       interpolated by the ``invalid_pixel_handler`` of `~ctapipe.calib.CameraCalibrator`.

    6. saturated pixels: the pixels whose maximum ADC sample (pedestal subtracted) is above
       ``saturation_threshold``, with more than ``saturation_width_threshold`` samples above
       ``saturation_width_level``, are flagged with `~ctapipe.containers.PixelStatus.SATURATED`
       in the pixel status of R1. Their charge is corrected after the image extraction by the
       `ImageSaturationCorrector`.

    The steps using the pedestal statistics (3 and the dead pixels of 5) are not applied
    if ``event.mon.tel[tel_id].r0`` is not filled.

    For the simulated events (``event.simulation`` filled), the R1 waveforms are given by the
    event source. They are only corrected for the PDE drop due to the NSB, and the
    ``bad_pixels`` are set to 0 and flagged (also in the true image).

    PDE drop correction (``pde_drop_factor``): the R1 waveforms are divided by this factor.
    The simulations use PDE files which include the drop for a given NSB level, while the
    observed data are corrected for it by the voltage drop correction: the factor must match
    the PDE file of the simulation (see ``mc_pde_correction_factors.json`` and the
    sst1mpipe_mc_config_low_nsb.json and sst1mpipe_mc_config_high_nsb.json configs).
    ``null`` (the default, and for the real telescopes 21 and 22) applies no correction.

    All the parameters can be set per telescope,
    e.g. ``"voltage_drop_correction": [["type", "*", "global"], ["id", 22, "none"]]``.
    """

    calibration_file = TelescopeParameter(
        trait=Path(exists=True, directory_ok=False, allow_none=True),
        default_value=None,
        allow_none=True,
        help=(
            "HDF5 file with the calibration parameters (dc_to_pe, calib_flag) from the"
            " dark runs. If None, the default calibration file of the telescope is used."
        ),
    ).tag(config=True)

    window_transmittance_file = TelescopeParameter(
        trait=Path(exists=True, directory_ok=False, allow_none=True),
        default_value=None,
        allow_none=True,
        help=(
            "File with the transmittance correction factor of the camera window of each pixel"
            " (pixel_id, correction_factor), measured in the lab: the charges are divided by it."
            " If None, the default file of the telescope (21, 22) is used."
        ),
    ).tag(config=True)

    voltage_drop_correction = TelescopeParameter(
        trait=CaselessStrEnum(VOLTAGE_DROP_CORRECTIONS),
        default_value="global",
        help=(
            "Correction of the voltage drop due to the NSB: 'global' uses the median"
            " of the pedestal variance of the camera, 'pixelwise' the variance of each pixel"
        ),
    ).tag(config=True)

    flag_bad_calibration_pixels = BoolTelescopeParameter(
        default_value=True,
        help=(
            "Set to 0 and flag the pixels with bad calibration parameters (calib_flag != 1)."
            " Otherwise they are calibrated with the mean dc_to_pe of the good pixels."
        ),
    ).tag(config=True)

    flag_dead_pixels = BoolTelescopeParameter(
        default_value=True,
        help="Set to 0 and flag the pixels whose pedestal std is below dead_pixel_std_threshold",
    ).tag(config=True)

    dead_pixel_std_threshold = FloatTelescopeParameter(
        default_value=2.5,
        help="Std of the ADC samples of the pedestal events (in ADC) below which a pixel is dead",
    ).tag(config=True)

    saturation_threshold = FloatTelescopeParameter(
        default_value=3000.0,
        help="ADC (pedestal subtracted) above which the maximum of a saturated pulse is",
    ).tag(config=True)

    saturation_width_level = FloatTelescopeParameter(
        default_value=2500.0,
        help="ADC (pedestal subtracted) level used to compute the width of the saturated pulses",
    ).tag(config=True)

    saturation_width_threshold = IntTelescopeParameter(
        default_value=5,
        help="Minimum number of samples above saturation_width_level of a saturated pulse",
    ).tag(config=True)

    bad_pixels = TelescopeParameter(
        trait=List(Int()),
        default_value=[("type", "*", [])],
        help=(
            "Ids of the pixels always set to 0 and flagged (e.g. broken pixels), for the"
            " observed and simulated events, e.g. [[\"type\", \"*\", []], [\"id\", 22, [3, 1201]]]"
        ),
    ).tag(config=True)

    pde_drop_factor = TelescopeParameter(
        trait=Float(allow_none=True),
        default_value=None,
        allow_none=True,
        help=(
            "The R1 waveforms are divided by this factor to correct the simulations for the"
            " PDE drop due to the NSB. It must match the PDE file of the simulation, see"
            " mc_pde_correction_factors.json. None: no correction (real telescopes)."
        ),
    ).tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        self._calibration = {}
        self._window_transmittance = {}
        self.n_bad_pixels = {}

    def calibration_file_path(self, tel_id):
        """Calibration file used for the telescope ``tel_id``"""
        try:
            path = self.calibration_file.tel[tel_id]
        except KeyError:
            path = None
        if path is not None:
            return path
        if tel_id not in DEFAULT_CALIBRATION_FILES:
            raise ValueError(f"No default calibration file for telescope {tel_id}, set calibration_file")
        return (CALIBRATION_DIR / DEFAULT_CALIBRATION_FILES[tel_id])

    def window_transmittance_file_path(self, tel_id):
        """Window transmittance file used for the telescope ``tel_id``"""
        try:
            path = self.window_transmittance_file.tel[tel_id]
        except KeyError:
            path = None
        if path is not None:
            return path
        if tel_id not in DEFAULT_WINDOW_FILES:
            raise ValueError(f"No default window transmittance file for telescope {tel_id}, set window_transmittance_file")
        return (WINDOW_DIR / DEFAULT_WINDOW_FILES[tel_id])

    def window_transmittance(self, tel_id):
        """Correction factor of the window transmittance of each pixel of the telescope ``tel_id``"""
        if tel_id not in self._window_transmittance:
            path = self.window_transmittance_file_path(tel_id)
            self._window_transmittance[tel_id] = read_window_transmittance(path)
            self.log.info("Telescope %d: window transmittance file %s", tel_id, path)
        return self._window_transmittance[tel_id]

    def calibration_parameters(self, tel_id):
        """
        dc_to_pe (with the mean of the good pixels for the bad pixels) and
        mask of the pixels with bad calibration parameters of the telescope ``tel_id``
        """
        if tel_id not in self._calibration:
            path = self.calibration_file_path(tel_id)
            parameters = pd.read_hdf(path)
            mask_bad = np.asarray(parameters['calib_flag'] != 1)
            dc_to_pe = np.array(parameters['dc_to_pe'], dtype=np.float64)
            dc_to_pe[mask_bad] = dc_to_pe[~mask_bad].mean()
            self._calibration[tel_id] = (dc_to_pe, mask_bad)
            self.log.info(
                "Telescope %d: calibration file %s, voltage drop correction: %s,"
                " flag bad calibration pixels: %s, flag dead pixels: %s (pedestal std < %.2f ADC)",
                tel_id, path, self.voltage_drop_correction.tel[tel_id],
                self.flag_bad_calibration_pixels.tel[tel_id], self.flag_dead_pixels.tel[tel_id],
                self.dead_pixel_std_threshold.tel[tel_id],
            )
        return self._calibration[tel_id]

    def voltage_drop(self, tel_id, pedestal_std):
        """Factor (scalar or per pixel) by which the charges are divided"""
        correction = self.voltage_drop_correction.tel[tel_id]
        if pedestal_std is None or correction == "none":
            return 1.0
        if correction == "global":
            return VAR_to_Idrop(np.median(pedestal_std**2), tel_id)
        return VAR_to_Idrop(pedestal_std**2, tel_id)

    def bad_pixel_mask(self, tel_id, pedestal_std):
        """Mask of the pixels set to 0 and flagged"""
        dc_to_pe, mask_bad_calibration = self.calibration_parameters(tel_id)
        mask_bad = np.zeros(dc_to_pe.shape, dtype=bool)

        if self.flag_bad_calibration_pixels.tel[tel_id]:
            mask_bad |= mask_bad_calibration
        if self.flag_dead_pixels.tel[tel_id] and pedestal_std is not None:
            mask_bad |= pedestal_std[0, :] < self.dead_pixel_std_threshold.tel[tel_id]
        mask_bad[self.static_bad_pixels(tel_id)] = True
        return mask_bad

    def static_bad_pixels(self, tel_id):
        """Ids of the ``bad_pixels`` of the telescope ``tel_id``"""
        try:
            return np.asarray(self.bad_pixels.tel[tel_id], dtype=int)
        except KeyError:
            return np.array([], dtype=int)

    def pde_drop(self, tel_id):
        """PDE drop factor by which the R1 waveforms are divided, None if not corrected"""
        try:
            return self.pde_drop_factor.tel[tel_id]
        except KeyError:
            return None

    def __call__(self, event, tel_id=None):
        """
        Calibrate the telescope ``tel_id`` of the event, or all its telescopes
        (with R0 data for the observed data, with R1 data for the simulations)
        """
        if event.simulation is not None:
            tel_ids = list(event.r1.tel.keys()) if tel_id is None else [tel_id]
            for tel in tel_ids:
                self._calibrate_simulated_telescope(event, tel)
            return event

        tel_ids = list(event.r0.tel.keys()) if tel_id is None else [tel_id]
        for tel in tel_ids:
            self._calibrate_telescope(event, tel)
        return event

    def _calibrate_simulated_telescope(self, event, tel_id):
        r1 = event.r1.tel[tel_id]
        pde_drop = self.pde_drop(tel_id)
        if pde_drop is not None:
            r1.waveform = r1.waveform / pde_drop

        bad_pixels = self.static_bad_pixels(tel_id)
        if len(bad_pixels) > 0:
            mask_bad = np.zeros(r1.waveform.shape[1], dtype=bool)
            mask_bad[bad_pixels] = True
            r1.waveform[:, mask_bad] = 0
            if event.simulation.tel[tel_id].true_image is not None:
                event.simulation.tel[tel_id].true_image[mask_bad] = 0
            self._flag_pixels(event, tel_id, mask_bad)
            self.n_bad_pixels[tel_id] = int(mask_bad.sum())

    @staticmethod
    def _flag_pixels(event, tel_id, mask_bad):
        """Flag the pixels, so that their charge is interpolated by the CameraCalibrator"""
        pixel_status = event.mon.tel[tel_id].pixel_status
        pixel_status.hardware_failing_pixels = mask_bad[np.newaxis]
        pixel_status.flatfield_failing_pixels = mask_bad[np.newaxis]
        pixel_status.pedestal_failing_pixels = mask_bad[np.newaxis]

    def _calibrate_telescope(self, event, tel_id):
        r0 = event.r0.tel[tel_id]
        dc_to_pe, _ = self.calibration_parameters(tel_id)
        pedestal_std = event.mon.tel[tel_id].r0.charge_std

        voltage_drop = np.asarray(self.voltage_drop(tel_id, pedestal_std))
        waveform = (r0.waveform - r0.pedestal[:, np.newaxis]) / dc_to_pe[:, np.newaxis]
        waveform /= voltage_drop[..., np.newaxis]
        waveform /= self.window_transmittance(tel_id)[:, np.newaxis]
        pde_drop = self.pde_drop(tel_id)
        if pde_drop is not None:
            waveform /= pde_drop

        # the R0 waveforms are kept: they are used afterwards by the R0 pedestal monitor
        mask_bad = self.bad_pixel_mask(tel_id, pedestal_std)
        waveform[:, mask_bad] = 0
        self.n_bad_pixels[tel_id] = int(mask_bad.sum())

        self._flag_pixels(event, tel_id, mask_bad)

        n_pixels = waveform.shape[1]
        pixel_status = np.full(n_pixels, PixelStatus.HIGH_GAIN_STORED, dtype=np.uint8)
        pixel_status[self.saturated_pixel_mask(tel_id, r0)] |= PixelStatus.SATURATED
        event.r1.tel[tel_id] = R1CameraContainer(
            event_type=event.trigger.event_type,
            event_time=event.trigger.tel[tel_id].time,
            waveform=waveform,
            selected_gain_channel=np.zeros(n_pixels, dtype=np.int8),
            pixel_status=pixel_status,
        )

    def saturated_pixel_mask(self, tel_id, r0):
        """Mask of the saturated pixels, from the ADC samples (pedestal subtracted)"""
        adc = r0.waveform[0] - r0.pedestal[:, np.newaxis]
        n_above_width_level = (adc > self.saturation_width_level.tel[tel_id]).sum(axis=1)
        return (
            (adc.max(axis=1) > self.saturation_threshold.tel[tel_id])
            & (n_above_width_level > self.saturation_width_threshold.tel[tel_id])
        )
