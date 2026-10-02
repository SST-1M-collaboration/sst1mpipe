import logging

import numpy as np
import pandas as pd
from importlib.resources import files

from ctapipe.containers import PixelStatus, R1CameraContainer
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import (
    BoolTelescopeParameter,
    CaselessStrEnum,
    FloatTelescopeParameter,
    Path,
    TelescopeParameter,
)

from sst1mpipe.utils import VAR_to_Idrop, get_tel_string



def get_default_window(telescope=None):
    """
    Provides default window transmissivity file,
    used in the case when it is not defined in the
    configuration file.

    Parameters
    ----------
    telescope: int
        Telescope number as in
        event.trigger.tels_with_trigger

    Returns
    -------
    window_corr: numpy.ndarray
    window_file: string

    """

    if (telescope == 21) or (telescope == 1):
        default_window_file_tel1 = 'corr_factor_1st_wdw.txt'
        logging.info('Window file used: ' + default_window_file_tel1)
        window_file = files('sst1mpipe.data').joinpath(default_window_file_tel1)
    elif (telescope == 22) or (telescope == 2):
        default_window_file_tel2 = 'corr_factor_2nd_wdw.txt'
        logging.info('Window file used: ' + default_window_file_tel2)
        window_file = files('sst1mpipe.data').joinpath(default_window_file_tel2)
    else:
        logging.error(f'Telescope {telescope} not known')
    window_corr = np.loadtxt(window_file, unpack=True, skiprows=1, usecols=1)
    return window_corr, window_file


def get_window_corr_factors(telescope=None, config=None):
    """
    Finds and reads the window transmissivity file.

    Parameters
    ----------
    telescope: int
        Telescope number as in
        event.trigger.tels_with_trigger
    config: dict

    Returns
    -------
    window_corr: numpy.ndarray
    window_file: string

    """

    if "window_transmittance" in config:
        if config["window_transmittance"]["tel_" + str(telescope).zfill(3)]:
            window_file = config["window_transmittance"]["tel_" + str(telescope).zfill(3)]
            window_corr = np.loadtxt(window_file, unpack=True, skiprows=1, usecols=1)
            logging.info("Window File for Tel %s: %s", telescope, window_file)
        else:
            logging.info("NO WINDOW TRANSMITTANCE FILE FOR TELESCOPE %s FOUND IN THE CFG FILE, DEFAULT WINDOW USED.", telescope)
            window_corr, window_file = get_default_window(telescope=telescope)
    else:
        logging.info("NO WINDOW TRANSMITTANCE FILE FOR TELESCOPE %s FOUND IN THE CFG FILE, DEFAULT WINDOW USED.", telescope)
        window_corr, window_file = get_default_window(telescope=telescope)

    return window_corr, window_file


def window_transmittance_correction(
        event, window_corr_factors=None,
        telescope=None,
        swapped_modules=None
        ):
    """
    Applies window transmittance correction
    on the integrated waveforms (images)

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer
    window_corr_factors: numpy.ndarray
    telescope: int
        Telescope number as in
        event.trigger.tels_with_trigger
    swapped_modules: list
        list of masks

    Returns
    -------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    """

    if swapped_modules is None:
        swapped_modules = []

    for mask_1,mask_2 in swapped_modules:

        # module 1
        window_corr_1 = window_corr_factors[mask_1]

        # module 2
        window_corr_2 = window_corr_factors[mask_2]

        window_corr_factors[mask_1] = window_corr_2
        window_corr_factors[mask_2] = window_corr_1

    image_corrected = event.dl1.tel[telescope].image / window_corr_factors
    event.dl1.tel[telescope].image = image_corrected.astype(np.float32)

    return event


def saturated_charge_correction(event):
    r"""
    Finds saturated waveforms and applies different peak integration on
    them, as the standard one does not perform well in such cases. This
    method integrates the peak above 20\% of the amplitude.
    Peak time for saturated events is also corrected as the middle of
    the integration window.

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    Returns
    -------
    saturated: bool
        True if the charges of saturated pixels were corrected

    """

    saturated_threshold = 3000
    width_level = 2500
    width_threshold = 5
    integration_level = 0.2

    telescope = event.trigger.tels_with_trigger[0]
    r0data = event.r0.tel[telescope]
    waveforms = (r0data.waveform[0].T - r0data.pedestal)

    # saturated pixels
    mask_saturated = np.max(waveforms, axis=0) > saturated_threshold
    saturated = False

    if sum(mask_saturated) > 0:

        image_new = event.dl1.tel[telescope].image
        peaktime_new = event.dl1.tel[telescope].peak_time

        # iterate over baseline subtracted waveforms and correct integration of those peaking above
        # saturation threshold and with larger width
        for k, w in enumerate(waveforms.T):

            if mask_saturated[k]:
                mask_width = w > width_level

                n_samples = mask_width.shape[0]

                max_adc = max(w)
                mask_integration = w >= integration_level * max_adc
                integration_start = np.arange(0, n_samples)[mask_integration][0]

                # This is needed to avoid secondary peaks (it looks for first drop below 20 percent after maximum)
                # We also select the first one, if there is a plato in the small peak in the wavefrom, which sometimes happen
                index_max = np.arange(0, n_samples)[w == max_adc][0]
                int_stop = np.arange(0, n_samples)[(w < integration_level * max_adc) & (np.arange(0, n_samples) > index_max)]
                # If it does not find where to stop it means the waveform extends over the readout window
                if len(int_stop) > 0:
                    integration_stop = int_stop[0]
                else:
                    integration_stop = n_samples-1

                width = sum(mask_width)
                peak_sample = width/2 + np.arange(0, n_samples)[mask_width][0]
                peak_time = peak_sample * 4

                if width > width_threshold:

                    # Peak integration correction
                    image_new[k] = sum(event.r1.tel[telescope].waveform[0, k][integration_start:integration_stop+1])
                    peaktime_new[k] = peak_time
                    saturated = True

        if saturated:
            event.dl1.tel[telescope].image = image_new
            event.dl1.tel[telescope].peak_time = peaktime_new

    return saturated


def correct_MC_for_PDE_drop(event, simtel_config_qe=None, pde_corr_factors=None):
    """
    Performs correction of the MC signal for the PDE (QE) drop.

    Parameters
    ----------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    simtel_config_qe: numpy.array
        Different PDEs (QEs) found in the simtel file header

    pde_corr_factors: dict
        Correction factors for different MC productions (it is NSB dependent)
        stored in ../data/mc_pde_correction_factors.json

    Returns
    -------
    event:
        sst1mpipe.io.containers.SST1MArrayEventContainer

    """

    for tel in event.r1.tel:

        try:
            stored_qe = pde_corr_factors['mc_correction_for_PDE'][get_tel_string(tel, mc=True)].keys()
            mask = [i in simtel_config_qe for i in stored_qe]
            VI = pde_corr_factors['mc_correction_for_PDE'][get_tel_string(tel, mc=True)][np.array(list(stored_qe))[mask][0]]
            event.r1.tel[tel].waveform /= VI
        except (KeyError, IndexError, TypeError, ZeroDivisionError):
            logging.error('PDE correction factors in the calibration file were not found in the simtel file header. Are you sure that you have listed the correct PDE factor in the mc_pde_correction_factors.json file?')
            exit()

    return event


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
    R0 -> R1 calibration of the SST-1M telescopes. Fills ``event.r1.tel[tel_id]``
    (`~ctapipe.containers.R1CameraContainer`, waveform in p.e. of shape
    (n_channels=1, n_pixels, n_samples)) from ``event.r0.tel[tel_id]``:

    1. subtraction of the pedestal computed by DigiCam (``r0.pedestal``)
    2. conversion from ADC to p.e. with the ``dc_to_pe`` of the calibration file
    3. voltage drop correction (``voltage_drop_correction``), from the std of the ADC
       samples of the pedestal events in ``event.mon.tel[tel_id].r0``
       (see `sst1mpipe.utils.monitoring_pedestals.R0PedestalMonitor`)
    4. bad pixels: pixels with bad calibration parameters (``flag_bad_calibration_pixels``)
       and dead pixels (``flag_dead_pixels``) are set to 0 in the R1 waveforms and
       flagged in ``event.mon.tel[tel_id].pixel_status``, so that their charge is
       interpolated by the ``invalid_pixel_handler`` of `~ctapipe.calib.CameraCalibrator`.

    The steps using the pedestal statistics (3 and the dead pixels of 4) are not applied
    if ``event.mon.tel[tel_id].r0`` is not filled. All the parameters can be set per
    telescope, e.g. ``"voltage_drop_correction": [["type", "*", "global"], ["id", 22, "none"]]``.
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

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray=subarray, config=config, parent=parent, **kwargs)
        self._calibration = {}
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
        return files('sst1mpipe.data').joinpath(DEFAULT_CALIBRATION_FILES[tel_id])

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

    def bad_pixels(self, tel_id, pedestal_std):
        """Mask of the pixels set to 0 and flagged"""
        dc_to_pe, mask_bad_calibration = self.calibration_parameters(tel_id)
        mask_bad = np.zeros(len(dc_to_pe), dtype=bool)
        if self.flag_bad_calibration_pixels.tel[tel_id]:
            mask_bad |= mask_bad_calibration
        if self.flag_dead_pixels.tel[tel_id] and pedestal_std is not None:
            mask_bad |= pedestal_std < self.dead_pixel_std_threshold.tel[tel_id]
        return mask_bad

    def __call__(self, event, tel_id=None):
        """
        Calibrate the telescope ``tel_id`` of the event, or all the telescopes with R0 data
        """
        tel_ids = list(event.r0.tel.keys()) if tel_id is None else [tel_id]
        for tel in tel_ids:
            self._calibrate_telescope(event, tel)
        return event

    def _calibrate_telescope(self, event, tel_id):
        r0 = event.r0.tel[tel_id]
        dc_to_pe, _ = self.calibration_parameters(tel_id)
        pedestal_std = event.mon.tel[tel_id].r0.charge_std

        voltage_drop = np.asarray(self.voltage_drop(tel_id, pedestal_std))
        waveform = (r0.waveform - r0.pedestal[:, np.newaxis]) / dc_to_pe[:, np.newaxis]
        waveform /= voltage_drop[..., np.newaxis]

        # the R0 waveforms are kept: they are used afterwards by the R0 pedestal monitor
        mask_bad = self.bad_pixels(tel_id, pedestal_std)
        waveform[:, mask_bad] = 0
        self.n_bad_pixels[tel_id] = int(mask_bad.sum())

        pixel_status = event.mon.tel[tel_id].pixel_status
        pixel_status.hardware_failing_pixels = mask_bad[np.newaxis]
        pixel_status.flatfield_failing_pixels = mask_bad[np.newaxis]
        pixel_status.pedestal_failing_pixels = mask_bad[np.newaxis]

        n_pixels = waveform.shape[1]
        event.r1.tel[tel_id] = R1CameraContainer(
            event_type=event.trigger.event_type,
            event_time=event.trigger.tel[tel_id].time,
            waveform=waveform,
            selected_gain_channel=np.zeros(n_pixels, dtype=np.int8),
            pixel_status=np.full(n_pixels, PixelStatus.HIGH_GAIN_STORED, dtype=np.uint8),
        )
