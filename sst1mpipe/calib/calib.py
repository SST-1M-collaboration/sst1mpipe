import logging

import astropy.units as u
import numpy as np
import pandas as pd

from ctapipe.containers import PixelStatus, R1CameraContainer
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import (
    BoolTelescopeParameter,
    Float,
    FloatTelescopeParameter,
    Int,
    IntTelescopeParameter,
    List,
    Path,
    TelescopeParameter,
)

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
    the maximum (or the end of the readout window), multiplied by the ``factor`` of the camera
    calibration coefficients (voltage drop correction), as the charges of the CameraCalibrator.
    Their peak time is the middle of the ADC samples (``event.r0``, pedestal subtracted) above
    ``peak_time_level``. The saturated pixels which are outliers (dead or bad pixels) keep the
    charge interpolated by the CameraCalibrator.

    ``integration_level`` or ``peak_time_level`` ``null``: the saturated pixels of the telescope are
    not corrected (simulated telescopes 1 and 2 in the MC configurations).
    """

    integration_level = TelescopeParameter(
        trait=Float(allow_none=True),
        default_value=0.2,
        allow_none=True,
        help=(
            "Fraction of the maximum defining the integration window of the saturated pulses."
            " None: no correction of the telescope"
        ),
    ).tag(config=True)

    peak_time_level = TelescopeParameter(
        trait=Float(allow_none=True),
        default_value=2500.0,
        allow_none=True,
        help=(
            "ADC (pedestal subtracted): the peak time of a saturated pulse is the middle"
            " of its R0 samples above this level. None: no correction of the telescope"
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
        """Pixels flagged as saturated in the pixel status of R1, without the outliers"""
        pixel_status = event.r1.tel[tel_id].pixel_status
        if pixel_status is None:
            return np.array([], dtype=int)
        saturated = (pixel_status & PixelStatus.SATURATED).astype(bool)
        outlier_mask = event.monitoring.tel[tel_id].camera.coefficients.outlier_mask
        if outlier_mask is not None:
            saturated &= ~np.atleast_2d(outlier_mask)[0]
        return np.flatnonzero(saturated)

    def _correct_telescope(self, event, tel_id):
        integration_level = self.integration_level.tel[tel_id]
        peak_time_level = self.peak_time_level.tel[tel_id]
        # the peak times are computed from the R0 waveforms: none for the R1 events of a stream
        if integration_level is None or peak_time_level is None or tel_id not in event.r0.tel:
            return False
        saturated = self.saturated_pixels(event, tel_id)
        if len(saturated) == 0:
            return False

        sample_time = (1 / self.subarray.tel[tel_id].camera.readout.sampling_rate).to_value(u.ns)
        factor = event.monitoring.tel[tel_id].camera.coefficients.factor
        factor = np.ones(len(event.dl1.tel[tel_id].image)) if factor is None else np.atleast_2d(factor)[0]

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
            dl1.image[pixel] = w[integration_start:integration_stop + 1].sum() * factor[pixel]

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


class R0R1Calibrator(TelescopeComponent):
    """
    R0 -> R1 calibration of the SST-1M telescopes: the corrections which do not need the pedestal
    events, so that it can be done without monitoring data (e.g. the ZMQ streams of the real-time
    analysis). For the observed data, fills ``event.r1.tel[tel_id]``
    (`~ctapipe.containers.R1CameraContainer`, waveform in p.e. of shape
    (n_channels=1, n_pixels, n_samples)) from ``event.r0.tel[tel_id]``:

    1. subtraction of the pedestal computed by DigiCam (``r0.pedestal``)
    2. conversion from ADC to p.e. with the ``dc_to_pe`` of the calibration file
    3. window transmittance correction: division by the correction factor of each pixel
       of ``window_transmittance_file``
    4. bad pixels: the pixels with bad calibration parameters (``flag_bad_calibration_pixels``)
       and the ``bad_pixels`` are added to ``event.monitoring.tel[tel_id].camera.coefficients.outlier_mask``,
       so that their charge is interpolated by the ``invalid_pixel_handler`` of
       `~ctapipe.calib.CameraCalibrator`
    5. saturated pixels: the pixels whose maximum ADC sample (pedestal subtracted) is above
       ``saturation_threshold``, with more than ``saturation_width_threshold`` samples above
       ``saturation_width_level``, are flagged with `~ctapipe.containers.PixelStatus.SATURATED`
       in the pixel status of R1. Their charge is corrected after the image extraction by the
       `ImageSaturationCorrector`.

    The corrections which need the pedestal events (voltage drop correction, dead pixels) are the
    camera calibration coefficients computed by sst1mpipe-monitoring
    (`~sst1mpipe.monitoring.CameraCoefficientsCalculator`), filled in
    ``event.monitoring.tel[tel_id].camera.coefficients`` by the monitoring source and applied by the
    CameraCalibrator. Without them (no monitoring source), they are not applied (warning).

    ``calibrate_r0`` false (by default for the simulated telescopes 1 and 2): the R1 waveforms
    given by the event source are kept, and the ``bad_pixels`` are added to the outliers (and set
    to 0 in the true image of the simulations).

    PDE drop correction (``pde_drop_factor``): the ``factor`` of the camera calibration coefficients
    is divided by it, i.e. the charges. The simulations use PDE files which include the drop for a
    given NSB level, while the observed data are corrected for it by the voltage drop correction:
    the factor must match the PDE file of the simulation (see ``mc_pde_correction_factors.json``
    and the sst1mpipe_mc_config_low_nsb.json and sst1mpipe_mc_config_high_nsb.json configs).
    ``null`` (the default, and for the real telescopes 21 and 22) applies no correction.

    All the parameters can be set per telescope,
    e.g. ``"bad_pixels": [["type", "*", []], ["id", 22, [3, 1201]]]``.
    """

    calibrate_r0 = BoolTelescopeParameter(
        default_value=[("type", "*", True), ("id", 1, False), ("id", 2, False)],
        help=(
            "Calibrate the R0 waveforms of the telescope (raw data). False: the R1 waveforms of the"
            " event source are kept (simulations, by default for the telescopes 1 and 2)"
        ),
    ).tag(config=True)

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

    flag_bad_calibration_pixels = BoolTelescopeParameter(
        default_value=True,
        help=(
            "Flag the pixels with bad calibration parameters (calib_flag != 1), whose charge is"
            " interpolated. Otherwise they are calibrated with the mean dc_to_pe of the good pixels."
        ),
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
            "Ids of the pixels always flagged (e.g. broken pixels), whose charge is interpolated, for"
            " the observed and simulated events, e.g. [[\"type\", \"*\", []], [\"id\", 22, [3, 1201]]]"
        ),
    ).tag(config=True)

    pde_drop_factor = TelescopeParameter(
        trait=Float(allow_none=True),
        default_value=None,
        allow_none=True,
        help=(
            "The charges are divided by this factor (factor of the camera calibration coefficients)"
            " to correct the simulations for the PDE drop due to the NSB. It must match the PDE file"
            " of the simulation, see mc_pde_correction_factors.json. None: no correction (real telescopes)."
        ),
    ).tag(config=True)

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        self._calibration = {}
        self._window_transmittance = {}
        self.n_bad_pixels = {}
        self._without_coefficients = set()

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
                "Telescope %d: calibration file %s, flag bad calibration pixels: %s",
                tel_id, path, self.flag_bad_calibration_pixels.tel[tel_id],
            )
        return self._calibration[tel_id]

    def bad_pixel_mask(self, tel_id, n_pixels):
        """Mask of the pixels flagged by the calibrator: bad calibration parameters (data) and ``bad_pixels``"""
        mask_bad = np.zeros(n_pixels, dtype=bool)
        if self.calibrate_r0.tel[tel_id] and self.flag_bad_calibration_pixels.tel[tel_id]:
            mask_bad |= self.calibration_parameters(tel_id)[1]
        mask_bad[self.static_bad_pixels(tel_id)] = True
        return mask_bad

    def static_bad_pixels(self, tel_id):
        """Ids of the ``bad_pixels`` of the telescope ``tel_id``"""
        try:
            return np.asarray(self.bad_pixels.tel[tel_id], dtype=int)
        except KeyError:
            return np.array([], dtype=int)

    def pde_drop(self, tel_id):
        """PDE drop factor by which the charges are divided, None if not corrected"""
        try:
            return self.pde_drop_factor.tel[tel_id]
        except KeyError:
            return None

    def __call__(self, event, tel_id=None):
        """
        Calibrate the telescope ``tel_id`` of the event, or all its telescopes: from their R0
        waveforms if ``calibrate_r0``, otherwise from the R1 waveforms of the event source
        """
        if tel_id is None:
            tel_ids = [tel for tel in event.r0.tel if self.calibrate_r0.tel[tel]]
            tel_ids += [tel for tel in event.r1.tel if not self.calibrate_r0.tel[tel]]
        else:
            tel_ids = [tel_id]
        for tel in tel_ids:
            if self.calibrate_r0.tel[tel]:
                self._calibrate_telescope(event, tel)
            n_pixels = event.r1.tel[tel].waveform.shape[1]
            mask_bad = self.bad_pixel_mask(tel, n_pixels)
            self.n_bad_pixels[tel] = int(mask_bad.sum())
            self._update_coefficients(event, tel, mask_bad)
            if event.simulation is not None and event.simulation.tel[tel].true_image is not None:
                event.simulation.tel[tel].true_image[mask_bad] = 0
        return event

    def _update_coefficients(self, event, tel_id, mask_bad):
        """
        Add the bad pixels to the outliers and the PDE drop to the factor of the camera calibration
        coefficients of the event (new arrays: the coefficients of the monitoring are not modified)
        """
        coefficients = event.monitoring.tel[tel_id].camera.coefficients
        if coefficients.outlier_mask is None:
            coefficients.outlier_mask = mask_bad[np.newaxis]
        else:
            coefficients.outlier_mask = np.asarray(coefficients.outlier_mask, dtype=bool) | mask_bad
        pde_drop = self.pde_drop(tel_id)
        if pde_drop is not None:
            factor = np.ones((1, len(mask_bad))) if coefficients.factor is None else coefficients.factor
            coefficients.factor = factor / pde_drop

    def _calibrate_telescope(self, event, tel_id):
        r0 = event.r0.tel[tel_id]
        if getattr(r0, "pedestal", None) is None:
            raise ValueError(
                f"Telescope {tel_id}: no DigiCam pedestal in the R0 data (simulation?),"
                " set R0R1Calibrator.calibrate_r0 to false for it"
            )
        if event.monitoring.tel[tel_id].camera.coefficients.factor is None and tel_id not in self._without_coefficients:
            self._without_coefficients.add(tel_id)
            self.log.warning(
                "Telescope %d: no camera calibration coefficients of the pedestal events (no monitoring"
                " source?): no voltage drop correction and no dead pixels", tel_id,
            )
        dc_to_pe, _ = self.calibration_parameters(tel_id)

        # the R0 waveforms are kept: they are used afterwards by the R0 pedestal monitor
        waveform = (r0.waveform - r0.pedestal[:, np.newaxis]) / dc_to_pe[:, np.newaxis]
        waveform /= self.window_transmittance(tel_id)[:, np.newaxis]

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
