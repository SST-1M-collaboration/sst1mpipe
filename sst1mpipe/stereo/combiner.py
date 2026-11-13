"""
Combination of the events of the SST-1M telescopes into stereo events.

The telescopes of SST-1M have no array trigger: each telescope writes its own events, and the
events of a shower seen by several telescopes have to be found afterwards. A `StereoEventCombiner`
does it in two steps:

1. `~StereoEventCombiner.match`: the events of the telescopes which are coincident, from their
   trigger tables (``/dl1/event/telescope/trigger`` of the DL1 files). Each child class implements a
   matching method (`~StereoEventCombiner.match_pair`):

   - `WhiteRabbitClosestCombiner`: closest White Rabbit trigger times (ns), within ``max_time_diff``
   - `SlidingWindowCombiner`: trigger times within a window, after a search of the offset between
     the clocks of the telescopes (data without White Rabbit timestamps)
   - `SWATEventIdCombiner`: same array event id given by SWAT

2. `~StereoEventCombiner.combine`: one event (`~sst1mpipe.io.containers.SST1MArrayEventContainer`)
   from the coincident events of the telescopes, one method by part of the event
   (``combine_index``, ``combine_trigger``, ``combine_r0``, ...), which a child class can override.

The events are matched to those of a reference telescope (``reference_tel_id``): a stereo event
is the event of the reference telescope and the events of the other telescopes coincident with it.
A coincidence of other telescopes without the reference telescope is not a stereo event.
Only the events of type SUBARRAY (triggered by the camera) are matched.
"""
from abc import abstractmethod

import astropy.units as u
import numpy as np
from ctapipe.containers import EventType
from ctapipe.core import Component
from ctapipe.core.traits import AstroQuantity, Int

from sst1mpipe.io.containers import SST1MArrayEventContainer

__all__ = [
    "SlidingWindowCombiner",
    "StereoEventCombiner",
    "SWATEventIdCombiner",
    "WhiteRabbitClosestCombiner",
    "nearest_index",
    "to_sst1m_event",
]


def _ns(quantity):
    """Integer number of ns of a time quantity"""
    return int(round(quantity.to_value(u.ns)))


def nearest_index(sorted_values, values):
    """Index of the nearest element of ``sorted_values`` (sorted, not empty) to each of ``values`` (the first one if two are as near)"""
    sorted_values, values = np.asarray(sorted_values), np.asarray(values)
    right = np.clip(np.searchsorted(sorted_values, values), 0, len(sorted_values) - 1)
    left = np.clip(right - 1, 0, len(sorted_values) - 1)
    left_closer = np.abs(values - sorted_values[left]) <= np.abs(sorted_values[right] - values)
    return np.where(left_closer, left, right)


def to_sst1m_event(event):
    """`SST1MArrayEventContainer` with the content of an event (e.g. read by the HDF5EventSource), unchanged"""
    if isinstance(event, SST1MArrayEventContainer):
        return event
    sst1m_event = SST1MArrayEventContainer()
    for name in event.fields:
        sst1m_event[name] = event[name]
    sst1m_event.meta = dict(event.meta)
    return sst1m_event


class StereoEventCombiner(Component):
    """
    Base class of the combination of the coincident events of the telescopes into stereo events:
    matching of the events (`match`, method of the child class) and combination of the matched
    events (`combine`).
    """

    reference_tel_id = Int(
        21,
        help=(
            "Reference telescope: the stereo events are its events with the coincident events of the"
            " other telescopes, and have its obs_id and event_id"
        ),
    ).tag(config=True)

    # ------------------------------------------------------------------ matching

    @abstractmethod
    def match_pair(self, reference, other):
        """
        Coincident events of the reference telescope and of another telescope.

        Parameters
        ----------
        reference, other : astropy.table.Table
            Events of the telescopes, sorted by trigger time: columns ``obs_id``, ``event_id`` and
            ``time_ns`` (trigger time, integer ns).

        Returns
        -------
        reference_index, other_index : numpy.ndarray
            Rows of the coincident events in each table, each row at most once.
        """

    def match(self, tables):
        """
        Stereo events: the coincident events of the telescopes.

        Parameters
        ----------
        tables : dict[int, astropy.table.Table]
            Events of each telescope (trigger tables): columns ``obs_id``, ``event_id``,
            ``event_type`` and ``time_ns`` (trigger time, integer ns), any order.

        Returns
        -------
        dict[int, numpy.ndarray]
            By telescope, the row in its table of the event of each stereo event, -1 for a
            telescope without event. The reference telescope has an event in each stereo event,
            and at least one other telescope.
        """
        if self.reference_tel_id not in tables:
            raise ValueError(f"No events of the reference telescope {self.reference_tel_id}")
        rows = {}
        for tel_id, table in tables.items():
            # only the events triggered by the camera, sorted by time
            selected = np.flatnonzero(np.asarray(table["event_type"]) == EventType.SUBARRAY.value)
            rows[tel_id] = selected[np.argsort(np.asarray(table["time_ns"])[selected], kind="stable")]

        reference_rows = rows[self.reference_tel_id]
        reference = tables[self.reference_tel_id][reference_rows]
        groups = {self.reference_tel_id: np.asarray(reference_rows)}
        for tel_id, table in tables.items():
            if tel_id == self.reference_tel_id:
                continue
            matched = np.full(len(reference_rows), -1)
            if len(reference_rows) > 0 and len(rows[tel_id]) > 0:
                reference_index, other_index = self.match_pair(reference, table[rows[tel_id]])
                matched[reference_index] = rows[tel_id][other_index]
            groups[tel_id] = matched
            self.log.info(
                "Telescope %d: %d events matched to the %d events of the reference telescope %d",
                tel_id, np.count_nonzero(matched >= 0), len(reference_rows), self.reference_tel_id,
            )
        stereo = np.zeros(len(reference_rows), dtype=bool)
        for tel_id, matched in groups.items():
            if tel_id != self.reference_tel_id:
                stereo |= matched >= 0
        return {tel_id: matched[stereo] for tel_id, matched in groups.items()}

    # ------------------------------------------------------------------ combination

    def combine(self, events):
        """
        Stereo event from the coincident events of the telescopes.

        Parameters
        ----------
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Event of each telescope (mono events, e.g. read by the HDF5EventSource).

        Returns
        -------
        SST1MArrayEventContainer
        """
        if self.reference_tel_id not in events:
            raise ValueError(f"A stereo event needs an event of the reference telescope {self.reference_tel_id}")
        stereo = SST1MArrayEventContainer()
        self.combine_index(stereo, events)
        self.combine_trigger(stereo, events)
        self.combine_r0(stereo, events)
        self.combine_r1(stereo, events)
        self.combine_dl0(stereo, events)
        self.combine_dl1(stereo, events)
        self.combine_dl2(stereo, events)
        self.combine_muon(stereo, events)
        self.combine_monitoring(stereo, events)
        # combine_array_pointing is not called: the merging of the array pointings is not defined yet
        if any(event.simulation is not None for event in events.values()):
            self.combine_simulation(stereo, events)
        return stereo

    def combine_index(self, stereo, events):
        """obs_id and event_id of the event of the reference telescope"""
        reference = events[self.reference_tel_id].index
        stereo.index.obs_id = reference.obs_id
        stereo.index.event_id = reference.event_id

    def combine_trigger(self, stereo, events):
        """
        Trigger of each telescope, all the telescopes in ``tels_with_trigger``, the earliest
        trigger time of the telescopes, and a SUBARRAY event
        """
        stereo.trigger.tels_with_trigger = np.array(sorted(events), dtype=np.int16)
        stereo.trigger.time = min(event.trigger.time for event in events.values())
        stereo.trigger.event_type = EventType.SUBARRAY
        for tel_id, event in events.items():
            stereo.trigger.tel[tel_id] = event.trigger.tel[tel_id]

    @staticmethod
    def _combine_telescope_map(stereo_level, events, level):
        """Telescope data of each event (``<level>.tel[tel_id]``) in the same level of the stereo event"""
        for tel_id, event in events.items():
            telescopes = event[level].tel
            if tel_id in telescopes:
                stereo_level.tel[tel_id] = telescopes[tel_id]

    def combine_r0(self, stereo, events):
        """R0 of each telescope"""
        self._combine_telescope_map(stereo.r0, events, "r0")

    def combine_r1(self, stereo, events):
        """R1 of each telescope"""
        self._combine_telescope_map(stereo.r1, events, "r1")

    def combine_dl0(self, stereo, events):
        """DL0 of each telescope"""
        self._combine_telescope_map(stereo.dl0, events, "dl0")

    def combine_dl1(self, stereo, events):
        """DL1 (images, parameters) of each telescope"""
        self._combine_telescope_map(stereo.dl1, events, "dl1")

    def combine_dl2(self, stereo, events):
        """DL2 of each telescope; the stereo DL2 is not filled (stereo reconstruction of another step)"""
        self._combine_telescope_map(stereo.dl2, events, "dl2")

    def combine_muon(self, stereo, events):
        """Muon parameters of each telescope"""
        self._combine_telescope_map(stereo.muon, events, "muon")

    def combine_monitoring(self, stereo, events):
        """Monitoring of each telescope (pointing, statistics...); the array pointing is not filled"""
        self._combine_telescope_map(stereo.monitoring, events, "monitoring")

    def combine_array_pointing(self, stereo, events):
        """
        Array pointing of the stereo event (``monitoring.pointing``): that of the events if they all
        have the same. Not called by `combine`.

        Raises
        ------
        NotImplementedError
            If the array pointings of the events differ: their merging is not defined yet.
        """
        pointings = [event.monitoring.pointing for event in events.values()]
        first = pointings[0]
        for pointing in pointings[1:]:
            for field in first.fields:
                a, b = u.Quantity(first[field]), u.Quantity(pointing[field])
                if not (np.all(a == b) or (np.all(np.isnan(a)) and np.all(np.isnan(b)))):
                    raise NotImplementedError("Combination of different array pointings is not implemented")
        stereo.monitoring.pointing = first

    def combine_simulation(self, stereo, events):
        """
        Raises
        ------
        NotImplementedError
            The simulated stereo events are already combined by sim_telarray.
        """
        raise NotImplementedError(
            "Combination of simulated events is not implemented: the stereo events of the simulations"
            " are combined by sim_telarray"
        )


class WhiteRabbitClosestCombiner(StereoEventCombiner):
    """
    Coincident events: closest White Rabbit trigger times (ns precision), within
    ``max_time_diff``. An event and its closest event of the other telescope are matched only if
    each one is the closest to the other: each event is matched at most once.
    """

    max_time_diff = AstroQuantity(
        default_value=10 * u.us,
        physical_type=u.physical.time,
        help="Maximum difference of the trigger times of coincident events (strict)",
    ).tag(config=True)

    def match_pair(self, reference, other):
        time_1 = np.asarray(reference["time_ns"], dtype=np.int64)
        time_2 = np.asarray(other["time_ns"], dtype=np.int64)
        index_1 = np.arange(len(time_1))
        index_2 = nearest_index(time_2, time_1)
        # mutual nearest events
        mutual = nearest_index(time_1, time_2)[index_2] == index_1
        close = np.abs(time_2[index_2] - time_1) < _ns(self.max_time_diff)
        keep = mutual & close
        return index_1[keep], index_2[keep]


class SlidingWindowCombiner(StereoEventCombiner):
    """
    Coincident events: trigger time of the other telescope within ``window_half_width`` of the
    trigger time of the reference telescope plus the offset between their clocks, for data without
    precise (White Rabbit) timestamps. The offset is ``time_offset`` if given, otherwise the offset
    giving the most coincident events, searched from ``offset_start`` to ``offset_stop`` in steps of
    ``offset_step``. An event with several events of the other telescope in its window is not
    matched (ambiguous), in both directions.
    """

    window_half_width = AstroQuantity(
        default_value=1 * u.ms, physical_type=u.physical.time,
        help="Half width of the coincidence window (strict)",
    ).tag(config=True)
    time_offset = AstroQuantity(
        default_value=None, allow_none=True, physical_type=u.physical.time,
        help="Offset of the clock of the other telescope (its time - reference time), searched if None",
    ).tag(config=True)
    offset_start = AstroQuantity(
        default_value=-10 * u.ms, physical_type=u.physical.time, help="First offset of the search",
    ).tag(config=True)
    offset_stop = AstroQuantity(
        default_value=10 * u.ms, physical_type=u.physical.time, help="End of the offsets of the search (excluded)",
    ).tag(config=True)
    offset_step = AstroQuantity(
        default_value=20 * u.us, physical_type=u.physical.time, help="Step of the offsets of the search",
    ).tag(config=True)

    def __init__(self, config=None, parent=None, **kwargs):
        super().__init__(config=config, parent=parent, **kwargs)
        #: offset (ns) used by the last call of match_pair
        self.last_time_offset = None

    @staticmethod
    def _window(sorted_times, times, half_width):
        """first index and number of ``sorted_times`` strictly within ``half_width`` of each of ``times``"""
        start = np.searchsorted(sorted_times, times - half_width, side="right")
        stop = np.searchsorted(sorted_times, times + half_width, side="left")
        return start, stop - start

    def n_coincidences(self, time_1, time_2, offset_ns):
        """Number of events of ``time_2`` within the window of the events of ``time_1`` + offset"""
        return int(self._window(time_2, time_1 + offset_ns, _ns(self.window_half_width))[1].sum())

    def find_time_offset(self, time_1, time_2):
        """Offset (ns) giving the most coincidences (the first one if several)"""
        offsets = np.arange(_ns(self.offset_start), _ns(self.offset_stop), _ns(self.offset_step), dtype=np.int64)
        if len(offsets) == 0:
            raise ValueError("SlidingWindowCombiner: no offset between offset_start and offset_stop")
        counts = [self.n_coincidences(time_1, time_2, offset) for offset in offsets]
        best = int(np.argmax(counts))
        self.log.info("Time offset found: %d ns (%d coincidences)", offsets[best], counts[best])
        return int(offsets[best])

    def match_pair(self, reference, other):
        time_1 = np.asarray(reference["time_ns"], dtype=np.int64)
        time_2 = np.asarray(other["time_ns"], dtype=np.int64)
        offset = _ns(self.time_offset) if self.time_offset is not None else self.find_time_offset(time_1, time_2)
        half_width = _ns(self.window_half_width)
        start_2, n_2 = self._window(time_2, time_1 + offset, half_width)
        _, n_1 = self._window(time_1, time_2 - offset, half_width)
        index_1 = np.flatnonzero(n_2 == 1)
        index_2 = start_2[index_1]
        # the event of the other telescope has a single event of the reference telescope in its window
        keep = n_1[index_2] == 1
        self.log.info(
            "%d ambiguous events (several events in the window)", np.count_nonzero(n_2 > 1) + np.count_nonzero(n_1 > 1),
        )
        self.last_time_offset = offset
        return index_1[keep], index_2[keep]


class SWATEventIdCombiner(StereoEventCombiner):
    """
    Coincident events: same array event id (``event_id``) given by SWAT. Only the event ids are
    used: events without the same id are not matched. An id found several times in a telescope
    is not matched.
    """

    def match_pair(self, reference, other):
        ids_1 = np.asarray(reference["event_id"])
        ids_2 = np.asarray(other["event_id"])
        _, index_1, index_2 = np.intersect1d(ids_1, ids_2, return_indices=True)
        values_1, counts_1 = np.unique(ids_1, return_counts=True)
        values_2, counts_2 = np.unique(ids_2, return_counts=True)
        repeated = np.union1d(values_1[counts_1 > 1], values_2[counts_2 > 1])
        keep = ~np.isin(ids_1[index_1], repeated)
        if np.any(~keep):
            self.log.warning("%d event ids found several times are not matched", np.count_nonzero(~keep))
        return index_1[keep], index_2[keep]
