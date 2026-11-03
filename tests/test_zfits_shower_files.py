from collections import Counter

import astropy.units as u
import numpy as np
import pytest
import tables
from astropy.time import Time
from ctapipe.containers import CoordinateFrameType, EventType, ObservingMode, PointingMode, SchedulingBlockType
from ctapipe.core import Provenance
from ctapipe.io import DataWriter, EventSource, read_table

from conftest import DARK_FILE, TRANSITION_FILES, run_process
from sst1mpipe.io import get_dl1_info
from sst1mpipe.io.containers import CameraEventType
import sst1mpipe.io.sst1m_event_source as sst1m_event_source
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.scripts.dark_run_mes_fitter import read_dark_run_events

N_PIXELS = 1296
N_SAMPLES = 50
N_TRIGGER_PATCHES = 432

# runs of the night of 2026-01-21 with Cherenkov events (camera trigger PATCH7)
# and pedestal events (internal trigger), with SWAT array event ids
FILES = {
    21: dict(
        path=TRANSITION_FILES[21], obs_id=202601210206,
        n_showers=84, n_pedestals=46, first_event_id=33759138, last_event_id=33759393,
        start="2026-01-21T18:49:45.441", stop="2026-01-21T18:49:45.894",
        date="2026-01-21T18:49:10", date_end="2026-01-21T18:49:10",
    ),
    22: dict(
        path=TRANSITION_FILES[22], obs_id=202601210585,
        n_showers=179, n_pedestals=67, first_event_id=43710705, last_event_id=43711231,
        start="2026-01-21T23:26:37.633", stop="2026-01-21T23:26:38.303",
        date="2026-01-21T23:26:02", date_end="2026-01-21T23:26:03",
    ),
}


@pytest.fixture(scope="module", params=sorted(FILES))
def zfits_file(request):
    tel_id = request.param
    return FILES[tel_id] | {"tel_id": tel_id}


@pytest.fixture(scope="module")
def events(zfits_file):
    """the fields of all the events of the file"""
    tel_id = zfits_file["tel_id"]
    summary = dict(
        obs_ids=set(), event_ids=[], times=[], tels_with_trigger=set(), event_types=Counter(),
        consistent_times=True, shapes=set(),
    )
    with SST1MEventSource(zfits_file["path"]) as source:
        summary["swat_event_ids_available"] = source.swat_event_ids_available
        for event in source:
            r0 = event.r0.tel[tel_id]
            summary["obs_ids"].add(event.index.obs_id)
            summary["event_ids"].append(event.index.event_id)
            summary["times"].append(event.trigger.time)
            summary["tels_with_trigger"].add(tuple(event.trigger.tels_with_trigger))
            summary["event_types"][(event.trigger.event_type, r0.event_type)] += 1
            summary["consistent_times"] &= bool(
                (r0.event_time == event.trigger.time) and (event.trigger.tel[tel_id].time == event.trigger.time)
            )
            summary["shapes"].add(tuple(
                (name, getattr(r0, name).shape, getattr(r0, name).dtype)
                for name in [
                    "waveform", "pedestal", "pixel_flags", "trigger_input_traces",
                    "trigger_output_patch7", "trigger_output_patch19", "trigger_output_muon",
                ]
            ))
    return summary


def test_zfits_file_is_read_by_sst1m_event_source(zfits_file):
    assert SST1MEventSource.is_compatible(zfits_file["path"])
    with EventSource(zfits_file["path"]) as source:
        assert isinstance(source, SST1MEventSource)
        assert not source.is_simulation
        assert source.obs_ids == [zfits_file["obs_id"]]
        # the telescope of the file only
        assert list(source.subarray.tel_ids) == [zfits_file["tel_id"]]


def test_event_ids(zfits_file, events):
    n_events = zfits_file["n_showers"] + zfits_file["n_pedestals"]

    assert events["obs_ids"] == {zfits_file["obs_id"]}
    assert events["tels_with_trigger"] == {(zfits_file["tel_id"],)}
    # SWAT array event ids: unique and increasing
    assert events["swat_event_ids_available"]
    event_ids = np.array(events["event_ids"])
    assert len(event_ids) == n_events
    assert np.all(np.diff(event_ids) > 0)
    assert event_ids[0] == zfits_file["first_event_id"]
    assert event_ids[-1] == zfits_file["last_event_id"]


def test_event_types(zfits_file, events):
    # Cherenkov events are triggered by the camera, pedestal events are internal triggers
    assert events["event_types"] == {
        (EventType.SUBARRAY, CameraEventType.PATCH7): zfits_file["n_showers"],
        (EventType.SKY_PEDESTAL, CameraEventType.INTERNAL): zfits_file["n_pedestals"],
    }


def test_event_times(zfits_file, events):
    times = Time(events["times"])
    assert times.scale == "tai"
    assert events["consistent_times"]
    assert np.all(np.diff(times.unix_tai) > 0)
    times.precision = 3
    assert times[0].isot == zfits_file["start"]
    assert times[-1].isot == zfits_file["stop"]


def test_r0_fields(events):
    assert events["shapes"] == {(
        ("waveform", (1, N_PIXELS, N_SAMPLES), np.dtype("int16")),
        ("pedestal", (N_PIXELS,), np.dtype("float64")),
        ("pixel_flags", (N_PIXELS,), np.dtype("uint16")),
        ("trigger_input_traces", (N_TRIGGER_PATCHES, N_SAMPLES), np.dtype("uint8")),
        ("trigger_output_patch7", (N_TRIGGER_PATCHES, N_SAMPLES), np.dtype("uint8")),
        ("trigger_output_patch19", (N_TRIGGER_PATCHES, N_SAMPLES), np.dtype("uint8")),
        ("trigger_output_muon", (N_TRIGGER_PATCHES, N_SAMPLES), np.dtype("uint8")),
    )}


def test_cherenkov_signal_in_waveforms(zfits_file):
    """
    a fraction of the Cherenkov events has a signal (maximum of the waveforms above the
    pedestal) larger than in any pedestal event: most of the triggered showers are small
    """
    tel_id = zfits_file["tel_id"]
    max_signal = {EventType.SUBARRAY: [], EventType.SKY_PEDESTAL: []}
    with SST1MEventSource(zfits_file["path"]) as source:
        for event in source:
            r0 = event.r0.tel[tel_id]
            signal = r0.waveform[0] - r0.pedestal[:, np.newaxis]
            max_signal[event.trigger.event_type].append(signal.max())

    showers = np.array(max_signal[EventType.SUBARRAY])
    assert np.mean(showers > max(max_signal[EventType.SKY_PEDESTAL])) > 0.1


def test_process_zfits_file(zfits_file, tmp_path):
    """sst1mpipe-process from R0 up to the DL1 parameters"""
    tel_id = zfits_file["tel_id"]
    output = tmp_path / "events.dl1.h5"
    # runs of the transition between two wobbles
    run_process(zfits_file["path"], output, "--ProcessorTool.allowed_sb_types=UNKNOWN")

    n_events = zfits_file["n_showers"] + zfits_file["n_pedestals"]
    trigger = read_table(output, "/dl1/event/subarray/trigger")
    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
    # the pedestal events are not written as events
    assert len(trigger) == len(parameters) == zfits_file["n_showers"]
    assert np.all(trigger["event_type"] == EventType.SUBARRAY.value)

    info = get_dl1_info(output)
    assert info["n_pedestal"][0] == zfits_file["n_pedestals"]
    assert info[f"n_triggered_tel{tel_id - 20}"][0] == n_events

    # images of Cherenkov events survive the cleaning
    intensity = parameters["camera_frame_hillas_intensity"]
    assert np.isfinite(intensity).sum() > 0
    assert np.all(intensity[np.isfinite(intensity)] > 0)


# ---------------------------------------------------------------------------
# scheduling and observation blocks of the run
# ---------------------------------------------------------------------------


def test_blocks_of_transition_run(zfits_file):
    """runs during the transition between two wobbles: scheduling block of UNKNOWN type"""
    obs_id = zfits_file["obs_id"]
    with SST1MEventSource(zfits_file["path"], max_events=1) as source:
        scheduling_block = source.scheduling_blocks[obs_id]
        observation_block = source.observation_blocks[obs_id]

    producer_id = f"SST1M-{zfits_file['tel_id']}"
    assert scheduling_block.sb_id == obs_id
    assert scheduling_block.sb_type == SchedulingBlockType.UNKNOWN
    assert scheduling_block.producer_id == producer_id
    assert scheduling_block.observing_mode == ObservingMode.UNKNOWN
    assert scheduling_block.pointing_mode == PointingMode.UNKNOWN

    assert observation_block.obs_id == obs_id
    assert observation_block.sb_id == obs_id
    assert observation_block.producer_id == producer_id
    assert observation_block.target == "Transition"
    assert observation_block.wobble == "NONE"
    assert np.isnan(observation_block.subarray_pointing_lon)
    assert observation_block.subarray_pointing_frame == CoordinateFrameType.UNKNOWN
    # DATE and DATEEND of the header of the file
    start, stop = Time(zfits_file["date"], scale="utc"), Time(zfits_file["date_end"], scale="utc")
    assert observation_block.actual_start_time == start
    assert u.isclose(observation_block.actual_duration, (stop - start).to(u.min))


def test_pointing_given_by_the_user():
    path = TRANSITION_FILES[21]
    with SST1MEventSource(path, max_events=1, pointing_ra=83.63, pointing_dec=22.01) as source:
        assert source.pointing_manual
        scheduling_block = source.scheduling_blocks[FILES[21]["obs_id"]]
        observation_block = source.observation_blocks[FILES[21]["obs_id"]]

    assert scheduling_block.pointing_mode == PointingMode.TRACK
    # the target of the file is kept
    assert scheduling_block.sb_type == SchedulingBlockType.UNKNOWN
    assert observation_block.target == "Transition"
    assert observation_block.subarray_pointing_lon.to_value(u.deg) == pytest.approx(83.63)
    assert observation_block.subarray_pointing_lat.to_value(u.deg) == pytest.approx(22.01)
    assert observation_block.subarray_pointing_frame == CoordinateFrameType.ICRS


def test_blocks_written_in_the_output(zfits_file, tmp_path):
    output = tmp_path / "events.dl1.h5"
    Provenance().start_activity("test_blocks_written_in_the_output")
    with SST1MEventSource(zfits_file["path"], max_events=1) as source:
        with DataWriter(source, output_path=output, write_dl1_parameters=True):
            pass

    scheduling_blocks = read_table(output, "/configuration/observation/scheduling_block")
    observation_blocks = read_table(output, "/configuration/observation/observation_block")
    assert list(scheduling_blocks["sb_id"]) == [zfits_file["obs_id"]]
    assert list(scheduling_blocks["sb_type"]) == [SchedulingBlockType.UNKNOWN.value]
    assert list(observation_blocks["obs_id"]) == [zfits_file["obs_id"]]
    assert list(observation_blocks["target"]) == ["Transition"]
    assert list(observation_blocks["wobble"]) == ["NONE"]
    # written in TAI by ctapipe
    written_start = observation_blocks["actual_start_time"][0]
    assert abs((written_start - Time(zfits_file["date"], scale="utc")).to_value(u.s)) < 1e-3


# ---------------------------------------------------------------------------
# runs which are not processed
# ---------------------------------------------------------------------------


def test_process_skips_transition_runs(zfits_file, tmp_path):
    """only the runs of the allowed scheduling block types (OBSERVATION by default) are processed"""
    output = tmp_path / "events.dl1.h5"
    tool = run_process(zfits_file["path"], output)

    n_events = zfits_file["n_showers"] + zfits_file["n_pedestals"]
    assert tool.n_skipped_events == n_events
    with tables.open_file(output) as h5:
        assert "/dl1/event" not in h5
    # the blocks and the production info are written
    scheduling_blocks = read_table(output, "/configuration/observation/scheduling_block")
    assert list(scheduling_blocks["sb_type"]) == [SchedulingBlockType.UNKNOWN.value]
    info = get_dl1_info(output)
    assert info["target"][0] == "Transition"
    assert info[f"n_triggered_tel{zfits_file['tel_id'] - 20}"][0] == 0


def test_allowed_sb_types_list(tmp_path):
    tool = run_process(
        TRANSITION_FILES[21], tmp_path / "events.dl1.h5", "--max-events=5",
        "--ProcessorTool.allowed_sb_types", "OBSERVATION", "--ProcessorTool.allowed_sb_types", "UNKNOWN",
    )

    assert tool.allowed_sb_types == [SchedulingBlockType.OBSERVATION, SchedulingBlockType.UNKNOWN]
    assert tool.n_skipped_events == 0


def test_dark_run_events(capsys):
    """the dark run script only uses the pedestal events of the dark runs"""
    events = list(read_dark_run_events([TRANSITION_FILES[21], DARK_FILE], max_events=20))

    assert len(events) == 20
    assert {event.index.obs_id for event in events} == {202601190007}
    assert {event.trigger.event_type for event in events} == {EventType.DARK_PEDESTAL}
    assert "is not a dark run" in capsys.readouterr().out


def test_dark_run_events_not_pedestal(monkeypatch, capsys):
    """the events of a dark run which are not dark pedestal events are not used"""
    # the run with Cherenkov events taken as a dark run
    monkeypatch.setattr(sst1m_event_source, "scheduling_block_type", lambda target: SchedulingBlockType.CALIBRATION)
    monkeypatch.setattr(sst1m_event_source, "pedestal_event_type", lambda target: EventType.DARK_PEDESTAL)
    path = TRANSITION_FILES[21]

    events = list(read_dark_run_events([path]))

    assert len(events) == FILES[21]["n_pedestals"]
    assert {event.trigger.event_type for event in events} == {EventType.DARK_PEDESTAL}
    assert f"{FILES[21]['n_showers']} events of the dark run" in capsys.readouterr().out
