import os.path

import pytest
from importlib.resources import files

from ctapipe.io import EventSource

from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.io.containers import CameraEventType, SST1MArrayEventContainer

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
CAMERA_EVENT_TYPE_1 = [CameraEventType.INTERNAL, CameraEventType.INTERNAL, CameraEventType.INTERNAL , CameraEventType.INTERNAL , CameraEventType.INTERNAL]

TEL_2_ID = 22

def test_test_tiles_exists():

    assert os.path.exists(FILE_TEL_1)
    assert os.path.exists(FILE_TEL_2)

def test_read_events():

    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=MAX_ITERATIONS)

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
        assert event.sst1m.r0.tel[TEL_1_ID].camera_event_type == CAMERA_EVENT_TYPE_1[i]
        i += 1
    assert i == MAX_ITERATIONS


def test_event_source_finds_sst1m_files():

    for path in (FILE_TEL_1, FILE_TEL_2):
        assert SST1MEventSource.is_compatible(path)
        with EventSource(input_url=path, max_events=1) as source:
            assert isinstance(source, SST1MEventSource)


def test_is_compatible_rejects_other_files():

    assert not SST1MEventSource.is_compatible(files('sst1mpipe.data').joinpath('sst1m_array.h5'))
    assert not SST1MEventSource.is_compatible(files('sst1mpipe.data').joinpath('sst1mpipe_data_config.json'))


@pytest.mark.parametrize("input_url", [
    [FILE_TEL_1, FILE_TEL_2],
    (FILE_TEL_1, FILE_TEL_2),
    [str(FILE_TEL_1), str(FILE_TEL_2)],
])
def test_input_url_list_of_files(input_url):

    source = SST1MEventSource(input_url=input_url, max_events=MAX_ITERATIONS)

    assert source.input_url == FILE_TEL_1
    assert source.filelist == [str(FILE_TEL_1), str(FILE_TEL_2)]

    # events are read starting with the first file
    event_ids = [event.sst1m.r0.event_id for event in source]
    assert event_ids == [FIRST_EVENT_ID_1 + i for i in range(MAX_ITERATIONS)]


def test_input_url_list_of_one_file():

    source = SST1MEventSource(input_url=[FILE_TEL_1])

    assert source.input_url == FILE_TEL_1
    assert source.filelist == [str(FILE_TEL_1)]


def test_files_are_read_one_after_the_other():

    n_events_file_1 = 11972
    source = SST1MEventSource(input_url=[FILE_TEL_1, FILE_TEL_2])

    for event in source:
        if event.sst1m.r0.tels_with_data[0] == TEL_2_ID:
            break

    assert event.count == n_events_file_1


def test_count_single_file():

    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=MAX_ITERATIONS)

    # count starts at 0 at each iteration over the source
    for _ in range(2):
        assert [event.count for event in source] == list(range(MAX_ITERATIONS))


@pytest.fixture
def fake_files(monkeypatch):
    """Replace the reading of a file by a few empty events, to test the loop over the files"""
    n_events = {str(FILE_TEL_1): 3, str(FILE_TEL_2): 2}

    def get_array_event(self, input_path):
        for _ in range(n_events[input_path]):
            event = SST1MArrayEventContainer()
            event.count = -1  # must be overwritten by the source
            event.meta["file"] = input_path
            yield event

    monkeypatch.setattr(SST1MEventSource, "get_array_event", get_array_event)
    return n_events


def test_count_continues_across_files(fake_files):

    source = SST1MEventSource(input_url=[FILE_TEL_1, FILE_TEL_2])
    events = [(event.count, event.meta["file"]) for event in source]

    assert events == [
        (0, str(FILE_TEL_1)), (1, str(FILE_TEL_1)), (2, str(FILE_TEL_1)),
        (3, str(FILE_TEL_2)), (4, str(FILE_TEL_2)),
    ]


def test_max_events_across_files(fake_files):

    source = SST1MEventSource(input_url=[FILE_TEL_1, FILE_TEL_2], max_events=4)
    events = [(event.count, event.meta["file"]) for event in source]

    assert [count for count, _ in events] == [0, 1, 2, 3]
    assert events[-1][1] == str(FILE_TEL_2)
