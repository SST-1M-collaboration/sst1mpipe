import logging
import os.path
import glob

from types import SimpleNamespace

import pytest

from ctapipe.core import Provenance
from ctapipe.io import EventSource, HDF5TableWriter, read_table

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time
from ctapipe.containers import CoordinateFrameType, ObservingMode, PointingMode, SchedulingBlockType

import sst1mpipe.io.sst1m_event_source as sst1m_event_source
from sst1mpipe.io.sst1m_event_source import (
    SST1MEventSource,
    file_has_swat_event_ids,
    file_start_and_duration,
    observing_mode,
    parse_file_name,
    parse_target_field,
    scheduling_block_type,
    tel_id_from_file_name,
)
from sst1mpipe.time import camera_clock_to_time
from sst1mpipe.io.containers import CameraEventType, DigicamConfigContainer, SST1MR0CameraContainer
from sst1mpipe.resources import DATA_CONFIG_FILE, SUBARRAY_FILE, TEST_DATA_DIR

# dark run of tel 22: pedestal events (internal triggers) only
DARK_FILE = (TEST_DATA_DIR / "zfits").joinpath('SST1M2_20260119_0007.fits.fz')
ZFITS_FILES = glob.glob(str(TEST_DATA_DIR / "zfits/*fits.fz"))

MAX_ITERATIONS = 5

DARK_TEL_ID = 22
OTHER_TEL_ID = 21
DARK_OBS_ID = 202601190007
DARK_FIRST_EVENT_ID = 28810122
DARK_FIRST_CAMERA_EVENT_NUMBER = 117266
DARK_SUM_WAVEFORM = [11926202, 11927313, 11927137, 11927213, 11927931]
DARK_SUM_PEDESTAL = [238504.25, 238498.625, 238495.125, 238505.5, 238511.625]
DARK_EVENT_TIME = [camera_clock_to_time(1768840365635653472 + i * 1_000_000) for i in range(MAX_ITERATIONS)]
DARK_CAMERA_EVENT_TYPE = [CameraEventType.INTERNAL] * MAX_ITERATIONS


def test_test_tiles_exists():

    assert os.path.exists(DARK_FILE)

def test_read_events():

    source = SST1MEventSource(input_url=DARK_FILE, max_events=MAX_ITERATIONS)

    i = 0
    for event in source:

        waveform = event.r0.tel[DARK_TEL_ID].waveform[0]
        baseline = event.r0.tel[DARK_TEL_ID].pedestal
        assert waveform.sum() == DARK_SUM_WAVEFORM[i]
        assert event.index.event_id == DARK_FIRST_EVENT_ID + i
        assert event.r0.tel[DARK_TEL_ID].camera_event_number == DARK_FIRST_CAMERA_EVENT_NUMBER + i
        assert baseline.sum() == DARK_SUM_PEDESTAL[i]
        assert event.r0.tel[DARK_TEL_ID].gps_time is None
        assert event.r0.tel[DARK_TEL_ID].event_time == DARK_EVENT_TIME[i]
        assert event.r0.tel[DARK_TEL_ID].event_type == DARK_CAMERA_EVENT_TYPE[i]
        i += 1
    assert i == MAX_ITERATIONS


def test_event_source_finds_sst1m_files():

    for path in ZFITS_FILES:

        assert SST1MEventSource.is_compatible(path)
        with EventSource(input_url=path, max_events=1) as source:
            assert isinstance(source, SST1MEventSource)


def test_is_compatible_rejects_other_files():

    assert not SST1MEventSource.is_compatible(SUBARRAY_FILE)
    assert not SST1MEventSource.is_compatible(DATA_CONFIG_FILE)


@pytest.mark.parametrize("input_url", [[DARK_FILE, DARK_FILE], (DARK_FILE,), [str(DARK_FILE)]])
def test_input_url_list_of_files_refused(input_url):
    """a single file is read: the processing scripts loop over the files"""
    with pytest.raises(TypeError, match="single file"):
        SST1MEventSource(input_url=input_url)


def test_count_single_file():

    source = SST1MEventSource(input_url=DARK_FILE, max_events=MAX_ITERATIONS)

    # count starts at 0 at each iteration over the source
    for _ in range(2):
        assert [event.count for event in source] == list(range(MAX_ITERATIONS))


@pytest.mark.parametrize("field, expected", [
    ("Crab_W1_83.63_22.01", ("Crab", "W1", 83.63, 22.01)),
    ("Crab,W2,83.63,22.01", ("Crab", "W2", 83.63, 22.01)),
    ("CrabW3_83.63_22.01", ("CrabW3", "W3", 83.63, 22.01)),
    ("Crab_83.63_22.01", ("Crab", "UNDEF", 83.63, 22.01)),
    ("MRK421_W1,166.994800,38.105300", ("MRK421", "W1", 166.9948, 38.1053)),
    ("Crab_W1_ra_dec", ("Crab", "W1", None, None)),
    ("Crab_W1_1_2_3", ("Crab", "W1", None, None)),
    ("dark", ("dark", None, None, None)),
    (None, (None, None, None, None)),
])
def test_parse_target_field(field, expected):

    assert parse_target_field(field) == expected


def test_trigger_and_no_pointing_for_dark_run():

    source = SST1MEventSource(input_url=DARK_FILE, max_events=MAX_ITERATIONS)

    assert source.target == "DARK"
    assert source.wobble is None
    assert source.pointing is None
    assert not source.pointing_manual
    observation_block = source.observation_blocks[DARK_OBS_ID]
    assert np.isnan(observation_block.subarray_pointing_lon)
    assert source.scheduling_blocks[DARK_OBS_ID].pointing_mode == PointingMode.UNKNOWN

    for i, event in enumerate(source):
        assert event.trigger.tels_with_trigger == [DARK_TEL_ID]
        assert event.trigger.time == DARK_EVENT_TIME[i]
        assert event.trigger.tel[DARK_TEL_ID].time == event.trigger.time
        assert np.isnan(event.monitoring.tel[DARK_TEL_ID].pointing.altitude)


@pytest.mark.parametrize("pointing_update_interval, tolerance", [(0, 1e-6 * u.arcsec), (3600, 1 * u.arcmin)])
def test_pointing_given_by_user(pointing_update_interval, tolerance):

    ra, dec = 83.633, 22.0145
    source = SST1MEventSource(
        input_url=DARK_FILE, max_events=MAX_ITERATIONS,
        pointing_ra=ra, pointing_dec=dec, pointing_update_interval=pointing_update_interval,
    )

    assert source.pointing_manual
    observation_block = source.observation_blocks[DARK_OBS_ID]
    assert observation_block.subarray_pointing_lon == ra * u.deg
    assert observation_block.subarray_pointing_lat == dec * u.deg
    assert observation_block.subarray_pointing_frame == CoordinateFrameType.ICRS
    assert source.scheduling_blocks[DARK_OBS_ID].pointing_mode == PointingMode.TRACK

    location = source.subarray.tel_coords.to_earth_location()[source.subarray.tel_index_array[DARK_TEL_ID]]
    target = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    for event in source:
        expected = target.transform_to(AltAz(obstime=event.trigger.time, location=location))
        pointing = event.monitoring.tel[DARK_TEL_ID].pointing
        filled = SkyCoord(az=pointing.azimuth, alt=pointing.altitude, frame=expected)
        assert filled.separation(expected) < tolerance
        assert u.isclose(event.monitoring.pointing.array_ra, ra * u.deg)
        assert u.isclose(event.monitoring.pointing.array_dec, dec * u.deg)


@pytest.mark.parametrize("target, expected", [
    ("Crab", SchedulingBlockType.OBSERVATION),
    ("Transition", SchedulingBlockType.UNKNOWN),
    ("TRANSITION", SchedulingBlockType.UNKNOWN),
    ("UNKNOWN", SchedulingBlockType.UNKNOWN),
    ("", SchedulingBlockType.UNKNOWN),
    (None, SchedulingBlockType.UNKNOWN),
    ("dark", SchedulingBlockType.CALIBRATION),
    ("DARK", SchedulingBlockType.CALIBRATION),
    ("drak", SchedulingBlockType.CALIBRATION),
    ("BIAS", SchedulingBlockType.CALIBRATION),
    ("WRtest", SchedulingBlockType.ENGINEERING),
])
def test_scheduling_block_type(target, expected):

    assert scheduling_block_type(target) == expected


@pytest.mark.parametrize("wobble, expected", [
    ("W1", ObservingMode.WOBBLE),
    ("W12", ObservingMode.WOBBLE),
    ("UNDEF", ObservingMode.UNKNOWN),
    (None, ObservingMode.UNKNOWN),
])
def test_observing_mode(wobble, expected):

    assert observing_mode(wobble) == expected


@pytest.mark.parametrize("file_name, expected", [
    ("SST1M1_20260121_0001.fits.fz", 21),
    ("/data/SST1M2_20251003_0123.fits.fz", 22),
    ("events.fits.fz", None),
])
def test_tel_id_from_file_name(file_name, expected):

    assert tel_id_from_file_name(file_name) == expected


def test_file_start_and_duration():

    start, duration = file_start_and_duration({"DATE": "2026-01-21T17:07:07", "DATEEND": "2026-01-21T17:07:19"})
    assert start == Time("2026-01-21T17:07:07", scale="utc")
    assert duration.to_value(u.s) == pytest.approx(12)
    assert file_start_and_duration({"DATE": "2026-01-21T17:07:07"}) == (None, None)


def test_blocks_of_dark_run():

    source = SST1MEventSource(input_url=DARK_FILE, max_events=1)
    scheduling_block = source.scheduling_blocks[DARK_OBS_ID]
    observation_block = source.observation_blocks[DARK_OBS_ID]

    assert scheduling_block.sb_id == DARK_OBS_ID
    assert scheduling_block.sb_type == SchedulingBlockType.CALIBRATION
    assert scheduling_block.producer_id == "SST1M-22"
    assert scheduling_block.observing_mode == ObservingMode.UNKNOWN
    assert scheduling_block.pointing_mode == PointingMode.UNKNOWN

    assert observation_block.obs_id == DARK_OBS_ID
    assert observation_block.sb_id == DARK_OBS_ID
    assert observation_block.producer_id == "SST1M-22"
    assert observation_block.target == "DARK"
    assert observation_block.wobble == "NONE"
    # DATE and DATEEND of the header of the file
    assert observation_block.actual_start_time == Time("2026-01-19T16:32:10", scale="utc")
    assert observation_block.actual_duration.to_value(u.s) == pytest.approx(34)


@pytest.mark.parametrize("file_name, expected", [
    ("SST1M1_20260121_0001.fits.fz", ("20260121", "0001")),
    ("/data/SST1M2_20251003_0123.fits.fz", ("20251003", "0123")),
    ("events.fits.fz", None),
])
def test_parse_file_name(file_name, expected):

    assert parse_file_name(file_name) == expected


def test_event_index():

    source = SST1MEventSource(input_url=DARK_FILE, max_events=MAX_ITERATIONS)

    assert source.run_id == DARK_OBS_ID
    assert list(source.observation_blocks) == [DARK_OBS_ID]
    assert source.observation_blocks[DARK_OBS_ID].obs_id == DARK_OBS_ID
    for i, event in enumerate(source):
        assert event.index.obs_id == DARK_OBS_ID
        assert event.index.event_id == DARK_FIRST_EVENT_ID + i


@pytest.mark.parametrize("input_file", ZFITS_FILES)
def test_swat_event_ids_in_files(input_file):

    assert file_has_swat_event_ids(input_file)
    assert SST1MEventSource(input_url=input_file, max_events=1).swat_event_ids_available


def test_swat_event_ids_used_as_event_id():

    source = SST1MEventSource(input_url=DARK_FILE, max_events=MAX_ITERATIONS)

    assert source.swat_event_ids_available
    for i, event in enumerate(source):
        # arrayEvtNum, not the camera event number
        assert event.index.event_id == DARK_FIRST_EVENT_ID + i
        assert event.index.event_id != event.r0.tel[DARK_TEL_ID].camera_event_number


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


def test_only_r0_trigger_and_pointing_are_filled():

    source = SST1MEventSource(input_url=DARK_FILE, max_events=MAX_ITERATIONS)
    n_pixels = source.subarray.tel[DARK_TEL_ID].camera.readout.n_pixels

    for event in source:
        assert list(event.r0.tel.keys()) == [DARK_TEL_ID]
        r0 = event.r0.tel[DARK_TEL_ID]
        assert isinstance(r0, SST1MR0CameraContainer)
        # one DigiCam channel
        assert r0.waveform.ndim == 3
        assert r0.pedestal.shape == (n_pixels, )
        assert "adc_samples" not in r0.fields
        assert "digicam_baseline" not in r0.fields
        assert r0.trigger_input_traces.shape == (432, r0.waveform.shape[-1])

        assert event.trigger.tels_with_trigger == [DARK_TEL_ID]
        assert event.trigger.time is not None

        assert len(event.r1.tel) == 0
        assert len(event.dl0.tel) == 0
        assert len(event.dl1.tel) == 0
        assert "sst1m" not in event.fields


def read_r0(**kwargs):
    source = SST1MEventSource(input_url=DARK_FILE, max_events=3, **kwargs)
    return [
        {key: getattr(event.r0.tel[DARK_TEL_ID], key).copy() for key in ("waveform", "pedestal", "pixel_flags")}
        for event in source
    ], source


# modules 8 and 9 of tel 21 wrongly connected during the run of the test file
SWAPPED_MODULES = [{"tel_id": DARK_TEL_ID, "start": "2026-01-19T00:00:00", "stop": "2026-01-20T00:00:00", "modules": [8, 9]}]


def test_no_swapped_modules_by_default():

    events, source = read_r0()

    assert source.swapped_pixels(DARK_TEL_ID, DARK_EVENT_TIME[0]) == []
    assert events[0]["waveform"].sum() == DARK_SUM_WAVEFORM[0]


def test_swap_modules():

    swapped, source = read_r0(swapped_modules=SWAPPED_MODULES)
    not_swapped, _ = read_r0()

    [(pixels_1, pixels_2)] = source.swapped_pixels(DARK_TEL_ID, DARK_EVENT_TIME[0])
    assert len(pixels_1) == len(pixels_2) == 12
    others = np.setdiff1d(np.arange(1296), np.concatenate([pixels_1, pixels_2]))
    # the telescope 22 is not affected
    assert source.swapped_pixels(OTHER_TEL_ID, DARK_EVENT_TIME[0]) == []

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

    source = SST1MEventSource(input_url=DARK_FILE, max_events=1, swapped_modules=SWAPPED_MODULES)
    day = 1 * u.day

    assert len(source.swapped_pixels(DARK_TEL_ID, DARK_EVENT_TIME[0])) == 1
    assert source.swapped_pixels(DARK_TEL_ID, DARK_EVENT_TIME[0] - day) == []
    assert source.swapped_pixels(DARK_TEL_ID, DARK_EVENT_TIME[0] + day) == []


def test_input_files_in_provenance():

    provenance = Provenance()
    provenance.start_activity("test_input_files_in_provenance")
    try:
        SST1MEventSource(input_url=DARK_FILE, max_events=1)
        inputs = provenance.current_activity.input
    finally:
        provenance.finish_activity()

    assert [entry["url"] for entry in inputs] == [str(DARK_FILE)]
    assert all(entry["role"] == "R0/Event" for entry in inputs)


def test_warning_no_pointing_in_file(caplog):

    with caplog.at_level(logging.WARNING):
        SST1MEventSource(input_url=DARK_FILE, max_events=1)

    assert "No pointing in the TARGET field ('DARK')" in caplog.text
    assert "not reconstructed" in caplog.text


def test_warning_pointing_given_by_user_and_in_file(monkeypatch, caplog):
    """the pointing of the file and the one given by the user are reported, with their separation"""
    header = {"TARGET": "Crab_W1_83.63_22.01", "DATE": "2026-01-21T17:07:07", "DATEEND": "2026-01-21T17:07:19"}
    monkeypatch.setattr(sst1m_event_source.fits, "getheader", lambda *args, **kwargs: header)

    with caplog.at_level(logging.WARNING):
        source = SST1MEventSource(input_url=DARK_FILE, max_events=1, pointing_ra=84.63, pointing_dec=22.01)

    assert source.pointing.ra.deg == pytest.approx(84.63)
    assert "Pointing given by the user (RA 84.6300 deg, Dec 22.0100 deg)" in caplog.text
    assert "pointing of the file (RA 83.6300 deg, Dec 22.0100 deg)" in caplog.text
    # 1 deg in RA at Dec 22 deg
    assert f"separation {np.cos(np.deg2rad(22.01)):.4f} deg" in caplog.text


def test_no_warning_pointing_in_file(monkeypatch, caplog):

    header = {"TARGET": "Crab_W1_83.63_22.01"}
    monkeypatch.setattr(sst1m_event_source.fits, "getheader", lambda *args, **kwargs: header)

    with caplog.at_level(logging.WARNING):
        source = SST1MEventSource(input_url=DARK_FILE, max_events=1)

    assert source.pointing.ra.deg == pytest.approx(83.63)
    assert source.scheduling_blocks[DARK_OBS_ID].sb_type == SchedulingBlockType.OBSERVATION
    assert source.scheduling_blocks[DARK_OBS_ID].observing_mode == ObservingMode.WOBBLE
    assert "pointing" not in caplog.text.lower()


@pytest.mark.parametrize("input_file, first_sn, digicam_time", [
    (DARK_FILE, 1120024, (626, 367976660)),
    (TEST_DATA_DIR / "zfits" / "SST1M1_20260120_1179.fits.fz", 2110003, (5122, 938643252)),
    (TEST_DATA_DIR / "zfits" / "SST1M2_20260121_0585.fits.fz", 1120024, (6865, 247889788)),
])
def test_digicam_config(input_file, first_sn, digicam_time):
    """configuration of the DigiCam boards, from the DigicamConfig table of the file"""
    source = SST1MEventSource(input_url=input_file, max_events=1)
    config = source.digicam_config

    assert isinstance(config, DigicamConfigContainer)
    # one entry per board slot, 0 for the empty slots
    for name in ["protocol_vers", "sn", "hv", "gateware_rev", "gateware_vers", "gateware_code",
                 "gateware_card_type", "firmware_rev", "firmware_vers", "firmware_code", "firmware_card_type"]:
        assert getattr(config, name).shape == (39,)
    assert config.sn.dtype == np.uint32
    boards = config.sn > 0
    assert boards.sum() == 34
    assert config.sn[boards][0] == first_sn
    assert np.all(config.protocol_vers[boards] == 1)
    assert set(config.firmware_rev[boards]) == {39, 46, 48}
    assert set(config.gateware_rev[boards]) == {23, 25}
    assert (config.digicam_time_sec, config.digicam_time_nanosec) == digicam_time
    assert (config.operation_id, config.operation_data) == (0, 0)


def test_digicam_config_can_be_written(tmp_path):

    config = SST1MEventSource(input_url=DARK_FILE, max_events=1).digicam_config
    with HDF5TableWriter(tmp_path / "config.h5") as writer:
        writer.write("digicam_config", config)

    table = read_table(tmp_path / "config.h5", "/digicam_config")
    np.testing.assert_array_equal(table["sn"][0], config.sn)
    assert table["digicam_time_sec"][0] == config.digicam_time_sec


def test_no_digicam_config(monkeypatch):

    class FileWithoutConfig:
        def __init__(self, path):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

    monkeypatch.setattr(sst1m_event_source, "File", FileWithoutConfig)
    assert sst1m_event_source.read_digicam_config("file.fits.fz") is None


@pytest.mark.parametrize("input_file, tel_id", [
    (TEST_DATA_DIR / "zfits" / "SST1M1_20260120_1179.fits.fz", 21),
    (TEST_DATA_DIR / "zfits" / "SST1M2_20260120_1102.fits.fz", 22),
])
def test_observation_run(input_file, tel_id, caplog):
    """observation of Mrk 421: blocks and pointing from the TARGET field of the file"""
    with caplog.at_level(logging.WARNING):
        source = SST1MEventSource(input_url=input_file, max_events=20)
    assert "pointing" not in caplog.text.lower()

    obs_id = source.run_id
    scheduling_block = source.scheduling_blocks[obs_id]
    observation_block = source.observation_blocks[obs_id]
    assert (source.target, source.wobble) == ("MRK421", "W1")
    assert not source.pointing_manual
    assert scheduling_block.sb_type == SchedulingBlockType.OBSERVATION
    assert scheduling_block.observing_mode == ObservingMode.WOBBLE
    assert scheduling_block.pointing_mode == PointingMode.TRACK
    assert scheduling_block.producer_id == f"SST1M-{tel_id}"
    assert (observation_block.target, observation_block.wobble) == ("MRK421", "W1")
    assert observation_block.subarray_pointing_lon.to_value(u.deg) == pytest.approx(166.9948)
    assert observation_block.subarray_pointing_lat.to_value(u.deg) == pytest.approx(38.1053)
    assert observation_block.subarray_pointing_frame == CoordinateFrameType.ICRS
    assert source.swat_event_ids_available

    # alt/az of the pointing of the telescope, recomputed every second
    location = source.subarray.tel_coords.to_earth_location()[source.subarray.tel_index_array[tel_id]]
    target = SkyCoord(ra=166.9948 * u.deg, dec=38.1053 * u.deg, frame="icrs")
    n_events = 0
    for event in source:
        expected = target.transform_to(AltAz(obstime=event.trigger.time, location=location))
        pointing = event.monitoring.tel[tel_id].pointing
        filled = SkyCoord(az=pointing.azimuth, alt=pointing.altitude, frame=expected)
        assert filled.separation(expected) < 1 * u.arcmin
        assert 74 < pointing.altitude.to_value(u.deg) < 75
        n_events += 1
    assert n_events == 20
