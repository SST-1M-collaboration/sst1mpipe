import os.path
from importlib.resources import files

from sst1mpipe.io.sst1m_event_source import SST1MEventSource

FILE_TEL_1 = files('sst1mpipe.resources.zfits').joinpath('SST1M1_20260121_0001.fits.fz')
FILE_TEL_2 = files('sst1mpipe.resources.zfits').joinpath('SST1M2_20260121_0001.fits.fz')

MAX_ITERATIONS = 5

TEL_1_ID = 21
FIRST_EVENT_ID_1 = 29901023
FIRST_CAMERA_EVENT_NUMBER_1 = 21079
SUM_WAVEFORM_1 = [22699322, 22695204, 22699302, 22700002, 22699551]
SUM_BASELINE_1 = [453938.875, 453933.25, 453942.0625, 453939.875, 453923.375]
LOCAL_CAMERA_CLOCK_1 = [1769015262383583536, 1769015262384583536, 1769015262385583536,
                        1769015262386583536, 1769015262387583536, ]

TEL_2_ID = 22

def test_test_tiles_exists():

    assert os.path.exists(FILE_TEL_1)
    assert os.path.exists(FILE_TEL_2)

def test_read_events():

    source = SST1MEventSource([str(FILE_TEL_1)], max_events=MAX_ITERATIONS)

    i = 0
    for event in source:

        waveform = event.sst1m.r0.tel[TEL_1_ID].adc_samples
        baseline = event.sst1m.r0.tel[TEL_1_ID].digicam_baseline
        assert waveform.sum() == SUM_WAVEFORM_1[i]
        assert event.sst1m.r0.event_id == FIRST_EVENT_ID_1 + i
        assert event.sst1m.r0.tel[TEL_1_ID].camera_event_number == FIRST_CAMERA_EVENT_NUMBER_1 + i
        assert baseline.sum() == SUM_BASELINE_1[i]
        assert event.sst1m.r0.tel[TEL_1_ID].gps_time == 0
        assert event.sst1m.r0.tel[TEL_1_ID].local_camera_clock == LOCAL_CAMERA_CLOCK_1[i]
        i += 1
    assert i == MAX_ITERATIONS
