import os.path

import logging
from types import SimpleNamespace

import pytest

from ctapipe.core import Provenance
from ctapipe.io import EventSource

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, SkyCoord
from ctapipe.containers import CoordinateFrameType, PointingMode

import sst1mpipe.io.sst1m_event_source as sst1m_event_source
from sst1mpipe.io.sst1m_event_source import (
    SST1MEventSource,
    file_has_swat_event_ids,
    parse_file_name,
    parse_target_field,
)
from sst1mpipe.time import camera_clock_to_time
from sst1mpipe.io import load_config
from sst1mpipe.io.containers import CameraEventType, SST1MArrayEventContainer, SST1MR0CameraContainer
from sst1mpipe.resources import DATA_CONFIG_FILE, SUBARRAY_FILE, TEST_DATA_DIR

FILE_TEL_1 = (TEST_DATA_DIR / "zfits").joinpath('SST1M1_20260121_0001.fits.fz')
FILE_TEL_2 = (TEST_DATA_DIR / "zfits").joinpath('SST1M2_20260121_0001.fits.fz')

MAX_ITERATIONS = 5

TEL_1_ID = 21
OBS_ID_1 = 202601210001
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

        waveform = event.r0.tel[TEL_1_ID].waveform[0]
        baseline = event.r0.tel[TEL_1_ID].pedestal
        assert waveform.sum() == SUM_WAVEFORM_1[i]
        assert event.index.event_id == FIRST_EVENT_ID_1 + i
        assert event.r0.tel[TEL_1_ID].camera_event_number == FIRST_CAMERA_EVENT_NUMBER_1 + i
        assert baseline.sum() == SUM_BASELINE_1[i]
        assert event.r0.tel[TEL_1_ID].gps_time == 0
        assert event.r0.tel[TEL_1_ID].event_time == LOCAL_CAMERA_CLOCK_1[i]
        assert event.r0.tel[TEL_1_ID].event_type == CAMERA_EVENT_TYPE_1[i]
        i += 1
    assert i == MAX_ITERATIONS


def test_event_source_finds_sst1m_files():

    for path in (FILE_TEL_1, FILE_TEL_2):
        assert SST1MEventSource.is_compatible(path)
        with EventSource(input_url=path, max_events=1) as source:
            assert isinstance(source, SST1MEventSource)


def test_is_compatible_rejects_other_files():

    assert not SST1MEventSource.is_compatible(SUBARRAY_FILE)
    assert not SST1MEventSource.is_compatible(DATA_CONFIG_FILE)


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
    event_ids = [event.index.event_id for event in source]
    assert event_ids == [FIRST_EVENT_ID_1 + i for i in range(MAX_ITERATIONS)]


def test_input_url_list_of_one_file():

    source = SST1MEventSource(input_url=[FILE_TEL_1])

    assert source.input_url == FILE_TEL_1
    assert source.filelist == [str(FILE_TEL_1)]


def test_files_are_read_one_after_the_other():

    n_events_file_1 = 11972
    source = SST1MEventSource(input_url=[FILE_TEL_1, FILE_TEL_2])

    for event in source:
        if event.trigger.tels_with_trigger[0] == TEL_2_ID:
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


@pytest.mark.parametrize("field, expected", [
    ("Crab_W1_83.63_22.01", ("Crab", "W1", 83.63, 22.01)),
    ("Crab,W2,83.63,22.01", ("Crab", "W2", 83.63, 22.01)),
    ("CrabW3_83.63_22.01", ("CrabW3", "W3", 83.63, 22.01)),
    ("Crab_83.63_22.01", ("Crab", "UNDEF", 83.63, 22.01)),
    ("Crab_W1_ra_dec", ("Crab", "W1", None, None)),
    ("Crab_W1_1_2_3", ("Crab", "W1", None, None)),
    ("dark", ("dark", None, None, None)),
    (None, (None, None, None, None)),
])
def test_parse_target_field(field, expected):

    assert parse_target_field(field) == expected


def test_camera_clock_to_time():

    time = camera_clock_to_time(LOCAL_CAMERA_CLOCK_1[1])

    assert time.scale == "tai"
    # ns precision is kept
    delta = time - camera_clock_to_time(LOCAL_CAMERA_CLOCK_1[0])
    assert np.isclose(delta.to_value(u.ns), LOCAL_CAMERA_CLOCK_1[1] - LOCAL_CAMERA_CLOCK_1[0], atol=1)


def test_trigger_and_no_pointing_for_dark_run():

    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=MAX_ITERATIONS)

    assert source.target == "dark"
    assert source.wobble is None
    assert source.pointing is None
    assert not source.pointing_manual
    observation_block = source.observation_blocks[OBS_ID_1]
    assert np.isnan(observation_block.subarray_pointing_lon)
    assert source.scheduling_blocks[OBS_ID_1].pointing_mode == PointingMode.UNKNOWN

    for i, event in enumerate(source):
        assert event.trigger.tels_with_trigger == [TEL_1_ID]
        assert event.trigger.time == camera_clock_to_time(LOCAL_CAMERA_CLOCK_1[i])
        assert event.trigger.tel[TEL_1_ID].time == event.trigger.time
        assert np.isnan(event.pointing.tel[TEL_1_ID].altitude)


@pytest.mark.parametrize("pointing_update_interval, tolerance", [(0, 1e-6 * u.arcsec), (3600, 1 * u.arcmin)])
def test_pointing_given_by_user(pointing_update_interval, tolerance):

    ra, dec = 83.633, 22.0145
    source = SST1MEventSource(
        input_url=FILE_TEL_1, max_events=MAX_ITERATIONS,
        pointing_ra=ra, pointing_dec=dec, pointing_update_interval=pointing_update_interval,
    )

    assert source.pointing_manual
    observation_block = source.observation_blocks[OBS_ID_1]
    assert observation_block.subarray_pointing_lon == ra * u.deg
    assert observation_block.subarray_pointing_lat == dec * u.deg
    assert observation_block.subarray_pointing_frame == CoordinateFrameType.ICRS
    assert source.scheduling_blocks[OBS_ID_1].pointing_mode == PointingMode.TRACK

    location = source.subarray.tel_coords.to_earth_location()[source.subarray.tel_index_array[TEL_1_ID]]
    target = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    for event in source:
        expected = target.transform_to(AltAz(obstime=event.trigger.time, location=location))
        pointing = event.pointing.tel[TEL_1_ID]
        filled = SkyCoord(az=pointing.azimuth, alt=pointing.altitude, frame=expected)
        assert filled.separation(expected) < tolerance
        assert u.isclose(event.pointing.array_ra, ra * u.deg)
        assert u.isclose(event.pointing.array_dec, dec * u.deg)


@pytest.mark.parametrize("file_name, expected", [
    ("SST1M1_20260121_0001.fits.fz", ("20260121", "0001")),
    ("/data/SST1M2_20251003_0123.fits.fz", ("20251003", "0123")),
    ("events.fits.fz", None),
])
def test_parse_file_name(file_name, expected):

    assert parse_file_name(file_name) == expected


def test_event_index():

    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=MAX_ITERATIONS)

    assert source.run_id == OBS_ID_1
    assert list(source.observation_blocks) == [OBS_ID_1]
    assert source.observation_blocks[OBS_ID_1].obs_id == OBS_ID_1
    for i, event in enumerate(source):
        assert event.index.obs_id == OBS_ID_1
        assert event.index.event_id == FIRST_EVENT_ID_1 + i


@pytest.mark.parametrize("input_file", [FILE_TEL_1, FILE_TEL_2])
def test_swat_event_ids_in_files(input_file):

    assert file_has_swat_event_ids(input_file)
    assert SST1MEventSource.check_swat_event_ids_available([input_file])
    assert SST1MEventSource.check_swat_event_ids_available(input_file)


def test_swat_event_ids_used_as_event_id():

    source = SST1MEventSource(input_url=[FILE_TEL_1, FILE_TEL_2], max_events=MAX_ITERATIONS)

    assert source.swat_event_ids_available
    for i, event in enumerate(source):
        # arrayEvtNum, not the camera event number
        assert event.index.event_id == FIRST_EVENT_ID_1 + i
        assert event.index.event_id != event.r0.tel[TEL_1_ID].camera_event_number


@pytest.fixture
def fake_array_event_numbers(monkeypatch):
    """Replace the zfits files by files with the given arrayEvtNum of their events"""
    array_event_numbers = {}

    class FakeFile:
        def __init__(self, path):
            self.Events = [SimpleNamespace(arrayEvtNum=n) for n in array_event_numbers[path]]

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

    monkeypatch.setattr(sst1m_event_source, "File", FakeFile)
    return array_event_numbers


@pytest.mark.parametrize("numbers, expected", [
    ([0, 0, 0, 0], False),  # not written by SWAT
    ([], False),  # empty file
    ([12], True),  # single event
    ([0, 12, 13], True),  # SWAT just restarted
    ([12, 0, 1], True),  # SWAT restarted during the file
    ([0] * 10 + [12], False),  # only the first events are read
])
def test_file_has_swat_event_ids(fake_array_event_numbers, numbers, expected):

    fake_array_event_numbers["file.fits.fz"] = numbers

    assert file_has_swat_event_ids("file.fits.fz") is expected


def test_file_has_swat_event_ids_n_events(fake_array_event_numbers):

    fake_array_event_numbers["file.fits.fz"] = [0, 0, 12]

    assert not file_has_swat_event_ids("file.fits.fz", n_events=2)
    assert file_has_swat_event_ids("file.fits.fz", n_events=3)


@pytest.mark.parametrize("numbers_1, numbers_2, expected", [
    ([1, 2], [3, 4], True),
    ([0, 0], [0, 0], False),
    ([1, 2], [0, 0], False),  # mixing ids would give inconsistent event ids
    ([0, 0], [3, 4], False),
])
def test_check_swat_event_ids_available_all_files(fake_array_event_numbers, numbers_1, numbers_2, expected):

    fake_array_event_numbers.update({"file_1.fits.fz": numbers_1, "file_2.fits.fz": numbers_2})

    available = SST1MEventSource.check_swat_event_ids_available(["file_1.fits.fz", "file_2.fits.fz"])

    assert available is expected


def test_check_swat_event_ids_warns_if_only_some_files(fake_array_event_numbers, caplog):

    fake_array_event_numbers.update({"file_1.fits.fz": [1, 2], "file_2.fits.fz": [0, 0]})

    with caplog.at_level(logging.WARNING):
        assert not SST1MEventSource.check_swat_event_ids_available(["file_1.fits.fz", "file_2.fits.fz"])

    assert "only in some of the files" in caplog.text


def test_check_swat_event_ids_no_file():

    assert not SST1MEventSource.check_swat_event_ids_available([])


def test_only_r0_trigger_and_pointing_are_filled():

    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=MAX_ITERATIONS)
    n_pixels = source.subarray.tel[TEL_1_ID].camera.readout.n_pixels

    for event in source:
        assert list(event.r0.tel.keys()) == [TEL_1_ID]
        r0 = event.r0.tel[TEL_1_ID]
        assert isinstance(r0, SST1MR0CameraContainer)
        # one DigiCam channel
        assert r0.waveform.shape == (1, n_pixels, r0.num_samples)
        assert r0.pedestal.shape == (n_pixels, )
        assert "adc_samples" not in r0.fields
        assert "digicam_baseline" not in r0.fields
        assert r0.trigger_input_traces.shape == (432, r0.num_samples)

        assert event.trigger.tels_with_trigger == [TEL_1_ID]
        assert event.trigger.time is not None

        assert len(event.r1.tel) == 0
        assert len(event.dl0.tel) == 0
        assert len(event.dl1.tel) == 0
        assert "sst1m" not in event.fields


def read_r0(**kwargs):
    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=3, **kwargs)
    return [
        {key: getattr(event.r0.tel[TEL_1_ID], key).copy() for key in ("waveform", "pedestal", "pixel_flags")}
        for event in source
    ], source


# modules 8 and 9 of tel 21 wrongly connected during the run of the test file
SWAPPED_MODULES = [{"tel_id": TEL_1_ID, "start": "2026-01-21T00:00:00", "stop": "2026-01-22T00:00:00", "modules": [8, 9]}]


def test_no_swapped_modules_by_default():

    events, source = read_r0()

    assert source.swapped_pixels(TEL_1_ID, LOCAL_CAMERA_CLOCK_1[0]) == []
    assert events[0]["waveform"].sum() == SUM_WAVEFORM_1[0]


def test_swapped_modules_of_data_config():

    config = load_config(DATA_CONFIG_FILE)
    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=1, config=config)

    assert [entry["modules"] for entry in source.swapped_modules] == [[59, 88], [8, 9]]
    # no wrongly connected modules of tel 21 at the date of the test file
    assert source.swapped_pixels(TEL_1_ID, LOCAL_CAMERA_CLOCK_1[0]) == []
    one_month = 30 * 24 * 3600 * 10**9
    assert len(source.swapped_pixels(TEL_1_ID, LOCAL_CAMERA_CLOCK_1[0] - one_month)) == 1
    assert source.swapped_pixels(22, LOCAL_CAMERA_CLOCK_1[0] - one_month) == []


def test_swap_modules():

    swapped, source = read_r0(swapped_modules=SWAPPED_MODULES)
    not_swapped, _ = read_r0()

    [(pixels_1, pixels_2)] = source.swapped_pixels(TEL_1_ID, LOCAL_CAMERA_CLOCK_1[0])
    assert len(pixels_1) == len(pixels_2) == 12
    others = np.setdiff1d(np.arange(1296), np.concatenate([pixels_1, pixels_2]))
    # the telescope 22 is not affected
    assert source.swapped_pixels(22, LOCAL_CAMERA_CLOCK_1[0]) == []

    for event, event_ref in zip(swapped, not_swapped, strict=True):
        # the pixels are swapped in all the R0 quantities with a value per pixel
        for key, value in event.items():
            # pixel axis first (waveform: (n_channels, n_pixels, n_samples))
            value, value_ref = (np.moveaxis(v, 1, 0) if key == "waveform" else v for v in (value, event_ref[key]))
            np.testing.assert_array_equal(value[pixels_1], value_ref[pixels_2])
            np.testing.assert_array_equal(value[pixels_2], value_ref[pixels_1])
            np.testing.assert_array_equal(value[others], value_ref[others])
        assert not np.array_equal(event["waveform"], event_ref["waveform"])


def test_swapped_modules_period():

    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=1, swapped_modules=SWAPPED_MODULES)
    day = 24 * 3600 * 10**9

    assert len(source.swapped_pixels(TEL_1_ID, LOCAL_CAMERA_CLOCK_1[0])) == 1
    assert source.swapped_pixels(TEL_1_ID, LOCAL_CAMERA_CLOCK_1[0] - day) == []
    assert source.swapped_pixels(TEL_1_ID, LOCAL_CAMERA_CLOCK_1[0] + day) == []


def test_input_files_in_provenance():

    provenance = Provenance()
    provenance.start_activity("test_input_files_in_provenance")
    try:
        SST1MEventSource(input_url=[FILE_TEL_1, FILE_TEL_2], max_events=1)
        inputs = provenance.current_activity.input
    finally:
        provenance.finish_activity()

    assert [entry["url"] for entry in inputs] == [str(FILE_TEL_1), str(FILE_TEL_2)]
    assert all(entry["role"] == "R0/Event" for entry in inputs)
