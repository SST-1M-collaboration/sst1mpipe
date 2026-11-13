"""
StereoEventCombiner components (sst1mpipe.stereo): matching of the coincident events of the
telescopes on made-up trigger tables, and combination of made-up events
"""
import astropy.units as u
import numpy as np
import pytest
from astropy.table import Table
from astropy.time import Time
from ctapipe.containers import ArrayEventContainer, EventType, ImageParametersContainer, SimulatedEventContainer
from traitlets.config import Config

from sst1mpipe.io.containers import SST1MArrayEventContainer
from sst1mpipe.stereo import (
    SlidingWindowCombiner,
    StereoEventCombiner,
    SWATEventIdCombiner,
    WhiteRabbitClosestCombiner,
    nearest_index,
    to_sst1m_event,
)

US = 1_000  # ns
MS = 1_000_000  # ns
T0 = 1_768_964_826 * 1_000_000_000  # ns


def triggers(tel_id, times_ns, event_ids=None, event_types=None, obs_id=None):
    """Trigger table of a telescope, as read by sst1mpipe-stereo"""
    n = len(times_ns)
    return Table({
        "obs_id": np.full(n, obs_id if obs_id is not None else 1000 + tel_id),
        "event_id": np.arange(n) if event_ids is None else np.asarray(event_ids),
        "tel_id": np.full(n, tel_id),
        "event_type": np.full(n, EventType.SUBARRAY.value) if event_types is None else np.asarray(event_types),
        "time_ns": T0 + np.asarray(times_ns, dtype=np.int64),
    })


def pairs(combiner, reference, other):
    """matched (reference row, other row) of two telescopes, with match()"""
    groups = combiner.match({21: reference, 22: other})
    return sorted(zip(groups[21].tolist(), groups[22].tolist(), strict=True))


def test_nearest_index():
    values = np.array([0, 10, 20, 30])
    np.testing.assert_array_equal(nearest_index(values, [-5, 4, 5, 6, 26, 100]), [0, 0, 0, 1, 3, 3])


def test_from_name_and_reference_config():
    config = Config({"StereoEventCombiner": {"reference_tel_id": 22}})
    for name, cls in [
        ("WhiteRabbitClosestCombiner", WhiteRabbitClosestCombiner),
        ("SlidingWindowCombiner", SlidingWindowCombiner),
        ("SWATEventIdCombiner", SWATEventIdCombiner),
    ]:
        combiner = StereoEventCombiner.from_name(name, config=config)
        assert isinstance(combiner, cls)
        # the configuration of the base class applies to all the combiners
        assert combiner.reference_tel_id == 22
    assert WhiteRabbitClosestCombiner().reference_tel_id == 21


# ---------------------------------------------------------------- WhiteRabbitClosestCombiner


def test_white_rabbit_closest():
    combiner = WhiteRabbitClosestCombiner()
    reference = triggers(21, [0, 100 * US, 200 * US, 300 * US])
    # 40 ns after 0; 9.99 us after 100 us; 10 us after 200 us (not < max_time_diff); none near 300 us
    other = triggers(22, [40, 100 * US + 9_990, 200 * US + 10 * US, 500 * US])
    assert pairs(combiner, reference, other) == [(0, 0), (1, 1)]
    combiner = WhiteRabbitClosestCombiner(max_time_diff=11 * u.us)
    assert pairs(combiner, reference, other) == [(0, 0), (1, 1), (2, 2)]


def test_white_rabbit_closest_one_to_one():
    """an event is matched at most once: with its closest event, if it is the closest to it"""
    combiner = WhiteRabbitClosestCombiner()
    # two events of tel 22 near the same event of tel 21: the closest one only
    assert pairs(combiner, triggers(21, [0]), triggers(22, [3 * US, 1 * US])) == [(0, 1)]
    # two events of tel 21 near the same event of tel 22
    assert pairs(combiner, triggers(21, [0, 5 * US]), triggers(22, [4 * US])) == [(1, 0)]


def test_white_rabbit_closest_not_sorted():
    """the rows are those of the tables given, in any order"""
    combiner = WhiteRabbitClosestCombiner()
    reference = triggers(21, [300 * US, 0, 100 * US])
    other = triggers(22, [100 * US + 50, 300 * US + 20])
    assert pairs(combiner, reference, other) == [(0, 1), (2, 0)]


def test_only_subarray_events_matched():
    combiner = WhiteRabbitClosestCombiner()
    pedestal = EventType.SKY_PEDESTAL.value
    subarray = EventType.SUBARRAY.value
    reference = triggers(21, [0, 100 * US], event_types=[pedestal, subarray])
    other = triggers(22, [10, 100 * US + 10], event_types=[subarray, pedestal])
    assert pairs(combiner, reference, other) == []
    other = triggers(22, [10, 100 * US + 10])
    assert pairs(combiner, reference, other) == [(1, 1)]


def test_no_events():
    combiner = WhiteRabbitClosestCombiner()
    assert pairs(combiner, triggers(21, []), triggers(22, [0])) == []
    assert pairs(combiner, triggers(21, [0]), triggers(22, [])) == []


# ---------------------------------------------------------------- SlidingWindowCombiner


def test_sliding_window_offset_search():
    """the offset of the clocks giving the most coincidences is found, then the events are matched"""
    rng = np.random.default_rng(1)
    times_1 = np.sort(rng.choice(np.arange(0, 10_000) * 3 * MS, size=500, replace=False))
    offset = 3_340 * US
    # 300 coincident events, and 200 events of tel 22 alone
    times_2 = np.sort(np.concatenate([
        times_1[:300] + offset + rng.integers(-100, 100, 300) * US // 10,
        rng.choice(np.arange(0, 10_000) * 3 * MS + MS, size=200, replace=False),
    ]))
    combiner = SlidingWindowCombiner(window_half_width=100 * u.us)
    matched = pairs(combiner, triggers(21, times_1), triggers(22, times_2))
    assert abs(combiner.last_time_offset - offset) <= 100 * US
    assert len(matched) == 300
    assert all(abs(times_2[j] - times_1[i] - offset) < 100 * US for i, j in matched)


def test_sliding_window_fixed_offset_and_ambiguous():
    combiner = SlidingWindowCombiner(time_offset=1 * u.ms, window_half_width=10 * u.us)
    # event 0: one candidate; event 1: two candidates (ambiguous); events 2 and 3: the same single
    # candidate of tel 22 (ambiguous for tel 22); event 4: candidate 10 us away (strict window)
    reference = triggers(21, [0, 100 * MS, 200 * MS, 200 * MS + 5 * US, 300 * MS])
    other = triggers(22, [
        1 * MS + 5 * US,
        101 * MS - 2 * US, 101 * MS + 2 * US,
        201 * MS + 3 * US,
        301 * MS + 10 * US,
    ])
    assert pairs(combiner, reference, other) == [(0, 0)]
    assert combiner.last_time_offset == 1 * MS


def test_sliding_window_no_offsets():
    combiner = SlidingWindowCombiner(offset_start=1 * u.ms, offset_stop=1 * u.ms)
    with pytest.raises(ValueError, match="no offset"):
        pairs(combiner, triggers(21, [0]), triggers(22, [0]))


# ---------------------------------------------------------------- SWATEventIdCombiner


def test_swat_event_ids():
    combiner = SWATEventIdCombiner()
    reference = triggers(21, [0, 100 * US, 200 * US], event_ids=[10, 11, 12])
    # same ids, any time; no fallback on the times: the event at the same time with another id is not matched
    other = triggers(22, [5 * MS, 100 * US, 200 * US + 10], event_ids=[12, 99, 10])
    assert pairs(combiner, reference, other) == [(0, 2), (2, 0)]


def test_swat_repeated_ids_not_matched():
    combiner = SWATEventIdCombiner()
    reference = triggers(21, [0, 1 * MS, 2 * MS], event_ids=[10, 11, 11])
    other = triggers(22, [0, 1 * MS, 3 * MS], event_ids=[10, 11, 12])
    assert pairs(combiner, reference, other) == [(0, 0)]


# ---------------------------------------------------------------- several telescopes


def test_match_three_telescopes():
    """events matched to the reference telescope; a coincidence without it is not a stereo event"""
    combiner = WhiteRabbitClosestCombiner()
    tables = {
        21: triggers(21, [0, 100 * US, 200 * US]),
        22: triggers(22, [10, 100 * US + 10, 300 * US]),
        23: triggers(23, [20, 250 * US, 300 * US + 5]),
    }
    groups = combiner.match(tables)
    rows = sorted(zip(groups[21].tolist(), groups[22].tolist(), groups[23].tolist(), strict=True))
    # 3 telescopes; 21 and 22; 21 alone and 22 + 23 without 21: not stereo events
    assert rows == [(0, 0, 0), (1, 1, -1)]


def test_reference_telescope():
    tables = {21: triggers(21, [0, 100 * US]), 22: triggers(22, [10])}
    groups = WhiteRabbitClosestCombiner(reference_tel_id=22).match(tables)
    assert groups[22].tolist() == [0] and groups[21].tolist() == [0]
    with pytest.raises(ValueError, match="reference telescope 23"):
        WhiteRabbitClosestCombiner(reference_tel_id=23).match(tables)


# ---------------------------------------------------------------- combination of the events


def mono_event(tel_id, time, obs_id, event_id, intensity):
    """Event of one telescope, as read from its DL1 file"""
    event = ArrayEventContainer()
    event.index.obs_id, event.index.event_id = obs_id, event_id
    event.trigger.time = time
    event.trigger.event_type = EventType.SUBARRAY
    event.trigger.tels_with_trigger = np.array([tel_id])
    event.trigger.tel[tel_id].time = time
    event.trigger.tel[tel_id].event_type = EventType.SUBARRAY
    event.dl1.tel[tel_id].parameters = ImageParametersContainer()
    event.dl1.tel[tel_id].parameters.hillas.intensity = intensity
    event.dl2.tel[tel_id]
    event.monitoring.tel[tel_id].pointing.altitude = 70 * u.deg
    event.monitoring.pointing.array_ra = 166.99 * u.deg
    event.monitoring.pointing.array_dec = 38.1 * u.deg
    return event


@pytest.fixture
def events():
    time = Time("2026-01-20T22:00:00", scale="tai")
    return {
        21: mono_event(21, time + 32 * u.ns, 202601201179, 25085157, 766.5),
        22: mono_event(22, time, 202601201102, 25085157, 900.5),
    }


def test_combine(events):
    stereo = WhiteRabbitClosestCombiner().combine(events)
    assert isinstance(stereo, SST1MArrayEventContainer)
    # ids of the reference telescope, time of the earliest telescope, SUBARRAY event
    assert (stereo.index.obs_id, stereo.index.event_id) == (202601201179, 25085157)
    assert stereo.trigger.time == events[22].trigger.time
    assert stereo.trigger.event_type == EventType.SUBARRAY
    assert stereo.trigger.tels_with_trigger.tolist() == [21, 22]
    for tel_id, event in events.items():
        assert stereo.trigger.tel[tel_id].time == event.trigger.time
        assert stereo.dl1.tel[tel_id] is event.dl1.tel[tel_id]
        assert stereo.dl2.tel[tel_id] is event.dl2.tel[tel_id]
        assert stereo.monitoring.tel[tel_id] is event.monitoring.tel[tel_id]
        assert tel_id not in stereo.r0.tel and tel_id not in stereo.r1.tel
    assert stereo.dl1.tel[22].parameters.hillas.intensity == 900.5
    # no stereo reconstruction, no array pointing (combine_array_pointing is not called)
    assert len(stereo.dl2.stereo.geometry) == 0
    assert np.isnan(stereo.monitoring.pointing.array_ra)


def test_combine_needs_the_reference_telescope(events):
    with pytest.raises(ValueError, match="reference telescope 23"):
        WhiteRabbitClosestCombiner(reference_tel_id=23).combine(events)


def test_combine_array_pointing(events):
    combiner = WhiteRabbitClosestCombiner()
    stereo = SST1MArrayEventContainer()
    combiner.combine_array_pointing(stereo, events)
    assert u.isclose(stereo.monitoring.pointing.array_ra, 166.99 * u.deg)
    events[22].monitoring.pointing.array_ra = 167 * u.deg
    with pytest.raises(NotImplementedError):
        combiner.combine_array_pointing(SST1MArrayEventContainer(), events)


def test_combine_simulation_not_implemented(events):
    events[21].simulation = SimulatedEventContainer()
    with pytest.raises(NotImplementedError, match="sim_telarray"):
        WhiteRabbitClosestCombiner().combine(events)


def test_combination_overridden_in_child_class(events):
    """a part of the combination is changed by overriding its method"""

    class LowestTelescopeIds(SWATEventIdCombiner):
        def combine_index(self, stereo, events):
            stereo.index = events[min(events)].index

    stereo = LowestTelescopeIds(reference_tel_id=22).combine(events)
    assert stereo.index.obs_id == 202601201179
    assert stereo.trigger.tels_with_trigger.tolist() == [21, 22]


def test_to_sst1m_event(events):
    event = to_sst1m_event(events[22])
    assert isinstance(event, SST1MArrayEventContainer)
    assert event.index is events[22].index and event.dl1 is events[22].dl1
    assert to_sst1m_event(event) is event
