from importlib.resources import files

import numpy as np
from ctapipe.calib import CameraCalibrator

from sst1mpipe.calib import Calibrator_R0_R1, saturated_charge_correction
from sst1mpipe.io import load_config
from sst1mpipe.io.sst1m_event_source import SST1MEventSource

FILE_TEL_1 = files('sst1mpipe.resources.zfits').joinpath('SST1M1_20260121_0001.fits.fz')
CONFIG = load_config(files('sst1mpipe.data').joinpath('sst1mpipe_data_config.json'), ismc=False)
TEL_ID = 21


def test_r0_r1_dl1_calibration():

    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=3)
    n_pixels = source.subarray.tel[TEL_ID].camera.readout.n_pixels
    calibrator_r0_r1 = Calibrator_R0_R1(config=CONFIG, telescope=TEL_ID)
    r1_dl1_calibrator = CameraCalibrator(subarray=source.subarray, config=CONFIG)

    for event in source:
        r0_waveform = event.r0.tel[TEL_ID].waveform.copy()

        calibrator_r0_r1.calibrate(event)
        r1 = event.r1.tel[TEL_ID]

        # ctapipe shape: (n_channels, n_pixels, n_samples)
        assert r1.waveform.shape == r0_waveform.shape == (1, n_pixels, r0_waveform.shape[-1])
        bad_pixels = event.mon.tel[TEL_ID].pixel_status.hardware_failing_pixels[0]
        assert bad_pixels.any()
        assert np.all(r1.waveform[:, bad_pixels] == 0)
        # the raw data are not modified by the calibration
        np.testing.assert_array_equal(event.r0.tel[TEL_ID].waveform, r0_waveform)

        r1_dl1_calibrator(event)
        assert event.dl1.tel[TEL_ID].image.shape == (n_pixels, )
        assert isinstance(saturated_charge_correction(event), bool)
