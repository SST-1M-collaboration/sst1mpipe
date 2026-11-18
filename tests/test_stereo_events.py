"""
Coincident (stereo) events of the two telescopes, in the observation of Mrk 421
of the test files: the run of tel 22 (35 events) is taken during the run of tel 21.
"""
import numpy as np
import pytest

from conftest import OBSERVATION_RUNS, read_zfits_events, run_process
from sst1mpipe.io import read_trigger_time
from sst1mpipe.io.containers import CameraEventType
from sst1mpipe.utils import get_trigger_time_ns

FILE_TEL_1, FILE_TEL_2 = OBSERVATION_RUNS[21]["path"], OBSERVATION_RUNS[22]["path"]
# events of tel 21 around the run of tel 22
EVENTS_TEL_1 = slice(4150, 4328)
# SWAT array event id -> trigger time of tel 22 - trigger time of tel 21 (ns)
COINCIDENT_EVENTS = {25085157: 32, 25085168: 48}
# maximum time difference of coincident events (stereo WhiteRabbitClosest method)
MAX_TIME_DIFF_NS = 10_000


def read_events(path, events=slice(None)):
    """SWAT event id, White Rabbit time (ns) and camera event type of the events of a zfits file"""
    events = read_zfits_events(path, events)
    return events["event_id"], events["time"], events["event_type"]


def match_closest_in_time(time_1, time_2, max_time_diff):
    """index of the closest event of the tel 1 for each event of tel 2, -1 if none within max_time_diff"""
    closest = np.argmin(np.abs(time_2[:, np.newaxis] - time_1[np.newaxis, :]), axis=1)
    time_diff = time_2 - time_1[closest]
    return np.where(np.abs(time_diff) < max_time_diff, closest, -1), time_diff


@pytest.fixture(scope="module")
def events_tel_1():
    return read_events(FILE_TEL_1, EVENTS_TEL_1)


@pytest.fixture(scope="module")
def events_tel_2():
    return read_events(FILE_TEL_2)


def test_runs_overlap(events_tel_1, events_tel_2):
    _, time_1, _ = events_tel_1
    _, time_2, _ = events_tel_2
    assert time_1[0] < time_2[0] and time_2[-1] < time_1[-1]


def test_coincident_events_swat_event_ids(events_tel_1, events_tel_2):
    """the coincident events have the same SWAT array event id"""
    ids_1, time_1, types_1 = events_tel_1
    ids_2, time_2, types_2 = events_tel_2

    common = np.intersect1d(ids_1, ids_2)
    assert set(common) == set(COINCIDENT_EVENTS)
    for event_id in common:
        index_1, index_2 = np.flatnonzero(ids_1 == event_id)[0], np.flatnonzero(ids_2 == event_id)[0]
        assert time_2[index_2] - time_1[index_1] == COINCIDENT_EVENTS[event_id]
        # triggered by the camera in both telescopes
        assert types_1[index_1] == types_2[index_2] == CameraEventType.PATCH7.value


def test_coincident_events_trigger_time(events_tel_1, events_tel_2):
    """the closest events in time (White Rabbit) are the events with the same SWAT event id"""
    ids_1, time_1, _ = events_tel_1
    ids_2, time_2, _ = events_tel_2

    matched, time_diff = match_closest_in_time(time_1, time_2, MAX_TIME_DIFF_NS)
    coincident = matched >= 0
    assert dict(zip(ids_2[coincident], time_diff[coincident], strict=True)) == COINCIDENT_EVENTS
    assert np.array_equal(ids_1[matched[coincident]], ids_2[coincident])
    # the other events of tel 22 are at least 100 us from an event of tel 21
    assert np.all(np.abs(time_diff[~coincident]) > 100_000)


def test_coincident_events_dl1(events_tel_1, tmp_path):
    """the trigger times of the DL1 file of tel 22 give the same coincident events"""
    output = tmp_path / "events.dl1.h5"
    run_process(FILE_TEL_2, output)
    trigger_time = read_trigger_time(output, "tel_022")

    ids_1, time_1, _ = events_tel_1
    matched, time_diff = match_closest_in_time(time_1, get_trigger_time_ns(trigger_time), MAX_TIME_DIFF_NS)
    coincident = matched >= 0
    assert dict(zip(trigger_time["event_id"][coincident], time_diff[coincident], strict=True)) == COINCIDENT_EVENTS
