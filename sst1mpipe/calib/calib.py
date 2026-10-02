import logging

import numpy as np
import pandas as pd
from importlib.resources import files

from ctapipe.containers import PixelStatus, R1CameraContainer

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


class Calibrator_R0_R1:

    """
    R0 -> R1 calibration

    Attributes
    ----------
        calibration_parameters: Pandas DataFrame
            Calibration parameters
        calibration_file: string
            File with calibration parameters
        dc_to_pe: numpy array
            ADC -> p.e. conversion factors
            (no gain drop)
        telescope: int
            Telescope number: 21/22
        config: dict
            Configuration
        mask_bad: numpy array
            Mask of bad pixels (True if bad)

    """

    def __init__(self, config=None, telescope=None):

        self.calibration_parameters = pd.DataFrame()
        self.calibration_file = ''
        self.dc_to_pe = np.array([])
        self.telescope = telescope
        self.config = config
        self.mask_bad = np.array([])
        self.pixels_removed = 0

        if self.config["telescope_calibration"]["bad_calib_px_interpolation"]:
            logging.info("Charge in pixels with bad calibration parameters are interpolated.")
            if self.config["telescope_calibration"]["dynamic_dead_px_interpolation"]:
                logging.info("Charge in dead pixels are interpolated (if pedestal info is provided).")
        else:
            logging.info("No pixel charge interpolation is applied. Gain in pixels with wrong calibration is calcualted as a global average of all gains.")

        # get calibration parameters during initialization
        self.get_calibration_parameters()
        self.get_dc_to_pe()


    def calibrate(self, event):

        """
        Runs the calibration. The voltage drop correction and the dead pixels
        identification use the statistics of the ADC samples of the pedestal events
        in event.mon.tel[tel].r0 (see `sst1mpipe.utils.monitoring_pedestals.R0PedestalMonitor`).
        They are not applied if this container is not filled.

        Parameters
        ----------
        event:
            sst1mpipe.io.containers.SST1MArrayEventContainer

        Returns
        -------
        event:
            sst1mpipe.io.containers.SST1MArrayEventContainer
        """

        r0data = event.r0.tel[self.telescope]
        baseline_subtracted = (r0data.waveform[0].T - r0data.pedestal)

        pedestal_std = event.mon.tel[self.telescope].r0.charge_std

        ## Apply (or not) pixel wise Voltage drop correction, VN: now also global
        ## TODO ?? TOTEST
        if pedestal_std is not None:
            if self.config['NsbCalibrator']['apply_pixelwise_Vdrop_correction']:
                VI = VAR_to_Idrop(pedestal_std**2, self.telescope)
            elif self.config['NsbCalibrator']['apply_global_Vdrop_correction']:
                VI = VAR_to_Idrop(np.median(pedestal_std**2), self.telescope)
            else:
                VI = 1.0
        else:
            VI = 1.0

        # as in ctapipe: waveform of shape (n_channels, n_pixels, n_samples),
        # with the single (high gain) channel of DigiCam
        n_pixels = baseline_subtracted.shape[1]
        event.r1.tel[self.telescope] = R1CameraContainer(
            event_type=event.trigger.event_type,
            event_time=event.trigger.tel[self.telescope].time,
            waveform=(baseline_subtracted / self.dc_to_pe / VI).T[np.newaxis],
            selected_gain_channel=np.zeros(n_pixels, dtype=np.int8),
            pixel_status=np.full(n_pixels, PixelStatus.HIGH_GAIN_STORED, dtype=np.uint8),
        )

        # This function removes bad pixels
        # - with not well determined dc_to_pe
        # - based on pedestal variation, i.e. dynamicaly
        # Charges in these pixels are then interpolated using method set in cfg: invalid_pixel_handler_type
        # Default is NeighborAverage, but can be turned off with 'null'
        if self.config["telescope_calibration"]["bad_calib_px_interpolation"]:
            event = self.remove_bad_pixels_calib(event)

        return event


    def get_calibration_parameters(self):

        """
        Finds and reads the calibration file.

        """

        if "telescope_calibration" in self.config:
            if self.config["telescope_calibration"]["tel_" + str(self.telescope).zfill(3)]:
                self.calib_file = self.config["telescope_calibration"]["tel_" + str(self.telescope).zfill(3)]
                self.calibration_parameters = pd.read_hdf(self.calib_file)
                logging.info("Calibration File for Tel %s: %s", self.telescope, self.calib_file)
            else:
                logging.info("NO CALIBRATION FILE FOR TELESCOPE %s FOUND IN THE CFG FILE, DEFAULT CALIBRATION USED.", self.telescope)
                self.get_default_calibration()
        else:
            logging.info("NO CALIBRATION FILE FOR TELESCOPE %s FOUND IN THE CFG FILE, DEFAULT CALIBRATION USED.", self.telescope)
            self.get_default_calibration()


    def get_default_calibration(self):
        """
        Provides default calibration file, used in
        the case when it is not defined in the
        configuration file.

        """

        # Calibration parameters averaged from all darks taken between
        # March 2023 and June 2024 (TEL1) and Sep 2023 and June 2024 (TEL2).
        # Based on TT's study, there is a relative difference in the dc_to_pe
        # factor between individual darks on the level of 5%, showing a
        # slowly decreasing trend. In the future, we may start producing
        # calibration files "per season" to mitigate the systematic uncertainty,
        # but per-night is not necessary. TT also confirms that dc_to_pe
        # does not depend on the level of DCR (the camera temperature).
        default_calib_file_tel1 = 'averaged_calib_param_v2_2023_2024_tel1.h5'
        default_calib_file_tel2 = 'averaged_calib_param_v2_2023_2024_tel2.h5'

        if (self.telescope == 21) or (self.telescope == 1):
            logging.info('Calib file used: ' + default_calib_file_tel1)
            self.calibration_file = files('sst1mpipe.data').joinpath(default_calib_file_tel1)
        elif (self.telescope == 22) or (self.telescope == 2):
            logging.info('Calib file used: ' + default_calib_file_tel2)
            self.calibration_file = files('sst1mpipe.data').joinpath(default_calib_file_tel2)
        else:
            logging.error(f'Telescope {self.telescope} not known')
        self.calibration_parameters = pd.read_hdf(self.calibration_file)


    def get_dc_to_pe(self):
        """
        Gets ADC -> p.e. conversion factor from the
        calibration file. Gain drop is not taken
        int account!

        """

        self.dc_to_pe = np.array(self.calibration_parameters['dc_to_pe'])
        # masking pixels with bad calibration parameters
        self.mask_bad = self.calibration_parameters['calib_flag'] != 1
        self.dc_to_pe[self.mask_bad] = self.dc_to_pe[~self.mask_bad].mean()


    def remove_bad_pixels_calib(self, event):
        """
        Fills bad pixel waveforms with zeros and
        flags them in proper containers. Charges in
        these pixels are then interpolated using method
        set in cfg: invalid_pixel_handler_type
        Default is NeighborAverage, but can be turned
        off with 'null'

        Parameters
        ----------
        event:
            sst1mpipe.io.containers.SST1MArrayEventContainer

        Returns
        -------
        event:
            sst1mpipe.io.containers.SST1MArrayEventContainer

        """

        # masking pixels with bad calibration parameters
        mask_bad_calib = self.mask_bad.astype(bool)

        # masking pixels with too low baseline std
        pedestal_std = event.mon.tel[self.telescope].r0.charge_std
        if pedestal_std is not None and self.config["telescope_calibration"]["dynamic_dead_px_interpolation"]:
            mask_bad_std = pedestal_std < 2.5
            mask_bad = mask_bad_calib + mask_bad_std
        else:
            mask_bad = mask_bad_calib
        self.pixels_removed = sum(mask_bad)

        # the R0 waveforms are kept: they are used afterwards by the R0 pedestal monitor
        event.r1.tel[self.telescope].waveform[:, mask_bad] = 0
        event.mon.tel[self.telescope].pixel_status['hardware_failing_pixels'] = np.array([mask_bad])
        event.mon.tel[self.telescope].pixel_status['flatfield_failing_pixels'] = np.array([mask_bad])
        event.mon.tel[self.telescope].pixel_status['pedestal_failing_pixels'] = np.array([mask_bad])

        return event
