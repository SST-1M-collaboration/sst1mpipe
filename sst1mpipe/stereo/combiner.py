"""
Combination of the events of the SST-1M telescopes into stereo events.

The SST-1M telescopes have no array trigger: each telescope triggers and writes its own events
(one raw data zfits file per run and per telescope, processed into one DL1 file each by
sst1mpipe-process), and the events of a shower seen by several telescopes
are found afterwards, from the trigger times or the array event ids of the events. A
`StereoEventCombiner` does it in two steps, used by the tool ``sst1mpipe-stereo``
(`~sst1mpipe.scripts.sst1mpipe_stereo_tool.StereoTool`):

1. `~StereoEventCombiner.match`: which events of the telescopes are coincident, from their
   trigger tables (``/dl1/event/telescope/trigger`` of the DL1 files). Each child class implements
   a matching method (`~StereoEventCombiner.match_pair`):

   - `WhiteRabbitClosestCombiner`: closest White Rabbit trigger times (ns), within ``max_time_diff``
   - `SlidingWindowCombiner`: trigger times within a window, after the search of the offset between
     the clocks of the telescopes (data without White Rabbit timestamps)
   - `SWATEventIdCombiner`: same array event id given by SWAT

2. `~StereoEventCombiner.combine`: one event (`~sst1mpipe.io.containers.SST1MArrayEventContainer`)
   from the coincident events of the telescopes, with one method per part of the event
   (`~StereoEventCombiner.combine_index`, `~StereoEventCombiner.combine_trigger`,
   `~StereoEventCombiner.combine_dl1`, ...), which a child class can override.

The events are matched to those of a reference telescope (``reference_tel_id``): a stereo event is
an event of the reference telescope with the events of the other telescopes coincident with it.
A coincidence of other telescopes without the reference telescope is not a stereo event. Only the
events of type SUBARRAY (triggered by the trigger of the camera: trigger patch7, trigger patch19 or
muon trigger) are matched, and each event is in at most one
stereo event.

Examples
--------
Combiner chosen by name (as by ``sst1mpipe-stereo --combiner``), matching of the events of two
telescopes and combination of the events of the first stereo event::

    from sst1mpipe.stereo import StereoEventCombiner

    combiner = StereoEventCombiner.from_name("WhiteRabbitClosestCombiner", reference_tel_id=21)
    groups = combiner.match({21: triggers_21, 22: triggers_22})
    # rows of the events of the first stereo event in the trigger tables
    row_21, row_22 = groups[21][0], groups[22][0]
    stereo_event = combiner.combine({21: event_21, 22: event_22})
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
    """
    Index of the nearest element of a sorted array to each of the given values.

    Parameters
    ----------
    sorted_values : array-like
        Values sorted in increasing order, not empty.
    values : array-like
        Values whose nearest element is searched.

    Returns
    -------
    numpy.ndarray
        For each of ``values``, the index of the nearest element of ``sorted_values``: the first
        one if two elements are as near.

    Examples
    --------
    >>> nearest_index([0, 10, 20], [4, 5, 6, 100])
    array([0, 0, 1, 2])
    """
    sorted_values, values = np.asarray(sorted_values), np.asarray(values)
    right = np.clip(np.searchsorted(sorted_values, values), 0, len(sorted_values) - 1)
    left = np.clip(right - 1, 0, len(sorted_values) - 1)
    left_closer = np.abs(values - sorted_values[left]) <= np.abs(sorted_values[right] - values)
    return np.where(left_closer, left, right)


def to_sst1m_event(event):
    """
    Event as a `~sst1mpipe.io.containers.SST1MArrayEventContainer`, with its content unchanged.

    Used for the mono events written by ``sst1mpipe-stereo --keep-mono``: the events read from the
    DL1 files by the ctapipe `~ctapipe.io.HDF5EventSource` are ctapipe
    `~ctapipe.containers.ArrayEventContainer`.

    Parameters
    ----------
    event : ctapipe.containers.ArrayEventContainer
        Event, e.g. read by the HDF5EventSource.

    Returns
    -------
    SST1MArrayEventContainer
        ``event`` itself if it is already a SST1MArrayEventContainer, otherwise a new container with
        the same fields (the same objects, not copies) and a copy of its ``meta``.
    """
    if isinstance(event, SST1MArrayEventContainer):
        return event
    sst1m_event = SST1MArrayEventContainer()
    for name in event.fields:
        sst1m_event[name] = event[name]
    sst1m_event.meta = dict(event.meta)
    return sst1m_event


class StereoEventCombiner(Component):
    """
    Base class of the combination of the coincident events of the telescopes into stereo events.

    A combiner matches the coincident events of the telescopes (`match`, from their trigger
    tables) and combines the matched events into a stereo event (`combine`). The child classes
    implement the matching of the events of two telescopes (`match_pair`), and can override the
    combination of any part of the event (``combine_<part>`` methods).

    The events of each telescope are matched to those of the reference telescope
    (``reference_tel_id``): a stereo event is an event of the reference telescope with the
    coincident events of the other telescopes, and has the obs_id and event_id of the event of
    the reference telescope. Only the events of type SUBARRAY are matched, and each event is in at
    most one stereo event.

    Parameters
    ----------
    config : traitlets.config.Config, optional
        Configuration of the combiner, e.g.
        ``Config({"WhiteRabbitClosestCombiner": {"max_time_diff": 5 * u.us}})``. The configuration
        of ``StereoEventCombiner`` (e.g. ``reference_tel_id``) applies to all the combiners.
    parent : ctapipe.core.Component or ctapipe.core.Tool, optional
        Parent of the combiner, whose configuration is used.
    **kwargs
        Values of the configuration options, e.g. ``reference_tel_id=22``.

    See Also
    --------
    WhiteRabbitClosestCombiner, SlidingWindowCombiner, SWATEventIdCombiner

    Warnings
    --------
    The combination of more than two telescopes is tested with three telescopes on made-up trigger
    tables only, and not with more than three telescopes: it may not be ready for such arrays.
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
        Coincident events of the reference telescope and of another telescope (method of the child class).

        Called by `match` for each telescope other than the reference telescope, with its events
        of type SUBARRAY only, sorted by trigger time. A child class implements its matching
        method here.

        Parameters
        ----------
        reference : astropy.table.Table
            Events of the reference telescope, sorted by trigger time, not empty: columns ``obs_id``,
            ``event_id``, ``tel_id``, ``event_type`` and ``time_ns`` (trigger time, integer ns).
        other : astropy.table.Table
            Events of the other telescope, sorted by trigger time, not empty, with the same columns.

        Returns
        -------
        reference_index : numpy.ndarray
            Rows of the coincident events in ``reference``.
        other_index : numpy.ndarray
            Rows of the coincident events in ``other``, in the same order: ``reference[reference_index[i]]``
            and ``other[other_index[i]]`` are coincident. Each row is at most once in each array.
        """

    def match(self, tables):
        """
        Stereo events: the coincident events of the telescopes.

        The events of type SUBARRAY of each telescope are matched to those of the reference
        telescope by `match_pair`. A stereo event is an event of the reference telescope with at
        least one coincident event of another telescope; with more than two telescopes, it has the
        coincident events of all the telescopes matched to the same event of the reference
        telescope.

        Parameters
        ----------
        tables : dict[int, astropy.table.Table]
            Events of each telescope by telescope id (e.g. the trigger tables of their DL1 files),
            in any order: columns ``obs_id``, ``event_id``, ``event_type`` (value of
            `~ctapipe.containers.EventType`) and ``time_ns`` (trigger time, integer ns).

        Returns
        -------
        dict[int, numpy.ndarray]
            By telescope id, the row in its table of its event in each stereo event (one element
            per stereo event, in the time order of the reference telescope), -1 if the telescope
            has no event in the stereo event. The reference telescope has an event in every stereo
            event.

        Raises
        ------
        ValueError
            If ``tables`` has no table of the reference telescope.
        """
        # TODO: test with more than three telescopes (and with three on real data)
        # TODO: match the coincidences of telescopes without an event of the reference telescope
        #       (matching of all the telescopes at once): they are mono events
        if self.reference_tel_id not in tables:
            raise ValueError(f"No events of the reference telescope {self.reference_tel_id}")
        rows = {}
        for tel_id, table in tables.items():
            # only the events triggered by the trigger of the camera (patch7, patch19 or muon), sorted by time
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

        Calls, in this order, `combine_index`, `combine_trigger`, `combine_r0`, `combine_r1`,
        `combine_dl0`, `combine_dl1`, `combine_dl2`, `combine_muon`, `combine_monitoring`, and
        `combine_simulation` if an event is simulated. `combine_array_pointing` is not called:
        the array pointing of the stereo event is not filled.

        The data of each telescope are those of its event: they are not copied, the stereo event
        refers to the same objects.

        Parameters
        ----------
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id: one event per telescope (mono events, e.g. read
            from the DL1 file of the telescope by the HDF5EventSource), with the event of the
            reference telescope.

        Returns
        -------
        SST1MArrayEventContainer
            Stereo event.

        Raises
        ------
        ValueError
            If ``events`` has no event of the reference telescope.
        NotImplementedError
            If an event is simulated (see `combine_simulation`).
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
        """
        Index of the stereo event (``stereo.index``): obs_id and event_id of the event of the reference telescope.

        All the ``combine_<part>`` methods have the same parameters, and fill their part of
        ``stereo`` in place.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id, with the event of the reference telescope.
        """
        reference = events[self.reference_tel_id].index
        stereo.index.obs_id = reference.obs_id
        stereo.index.event_id = reference.event_id

    def combine_trigger(self, stereo, events):
        """
        Trigger of the stereo event (``stereo.trigger``).

        - ``tels_with_trigger``: the telescopes of the events, sorted
        - ``time``: the earliest trigger time of the events
        - ``event_type``: SUBARRAY
        - ``tel[tel_id]``: the trigger of each telescope (time, event type, trigger pixels), from its event

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.
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
        """
        R0 of each telescope (``stereo.r0.tel[tel_id]``), from its event if it has R0 data.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.
        """
        self._combine_telescope_map(stereo.r0, events, "r0")

    def combine_r1(self, stereo, events):
        """
        R1 of each telescope (``stereo.r1.tel[tel_id]``), from its event if it has R1 data.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.
        """
        self._combine_telescope_map(stereo.r1, events, "r1")

    def combine_dl0(self, stereo, events):
        """
        DL0 of each telescope (``stereo.dl0.tel[tel_id]``), from its event if it has DL0 data.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.
        """
        self._combine_telescope_map(stereo.dl0, events, "dl0")

    def combine_dl1(self, stereo, events):
        """
        DL1 of each telescope (``stereo.dl1.tel[tel_id]``: image, peak time, image mask, parameters), from its event.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.
        """
        self._combine_telescope_map(stereo.dl1, events, "dl1")

    def combine_dl2(self, stereo, events):
        """
        DL2 of each telescope (``stereo.dl2.tel[tel_id]``), from its event.

        The stereo DL2 (``stereo.dl2.stereo``) is not filled: the stereo reconstruction (e.g. the
        geometry of the shower with the ctapipe `~ctapipe.reco.ShowerProcessor`) is a later step.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.
        """
        self._combine_telescope_map(stereo.dl2, events, "dl2")

    def combine_muon(self, stereo, events):
        """
        Muon parameters of each telescope (``stereo.muon.tel[tel_id]``), from its event.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.
        """
        self._combine_telescope_map(stereo.muon, events, "muon")

    def combine_monitoring(self, stereo, events):
        """
        Monitoring of each telescope (``stereo.monitoring.tel[tel_id]``: pointing, camera monitoring), from its event.

        The array pointing (``stereo.monitoring.pointing``) is not filled, see
        `combine_array_pointing`.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.
        """
        self._combine_telescope_map(stereo.monitoring, events, "monitoring")

    def combine_array_pointing(self, stereo, events):
        """
        Array pointing of the stereo event (``stereo.monitoring.pointing``). Not called by `combine`.

        The array pointing of the events (filled by the HDF5EventSource from the observation block
        of their run) is kept if all the events have the same. The combination of different array
        pointings is not defined yet.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event, filled in place.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.

        Raises
        ------
        NotImplementedError
            If the array pointings of the events differ (NaN values are equal).
        """
        # TODO: define the combination of different array pointings, then call this method in combine()
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
        Simulation of the stereo event: not implemented.

        Called by `combine` if an event is simulated. The stereo events of the simulations are
        already combined by sim_telarray (array trigger), with the same event_id in all the
        telescopes: they do not need to be matched.

        Parameters
        ----------
        stereo : SST1MArrayEventContainer
            Stereo event.
        events : dict[int, ctapipe.containers.ArrayEventContainer]
            Coincident events by telescope id.

        Raises
        ------
        NotImplementedError
            Always.
        """
        raise NotImplementedError(
            "Combination of simulated events is not implemented: the stereo events of the simulations"
            " are combined by sim_telarray"
        )


class WhiteRabbitClosestCombiner(StereoEventCombiner):
    """
    Coincident events: closest White Rabbit trigger times.

    For data with White Rabbit timestamps (ns precision, data taken after 25 September 2023). Each
    event of the reference telescope is matched to the event of the other telescope with the
    closest trigger time, if their trigger times differ by less than ``max_time_diff``, and if it
    is also the closest event of the reference telescope to that event (mutual closest events):
    each event is matched at most once.

    The trigger times of the coincident events of the SST-1M telescopes differ by a few tens of ns
    (e.g. 32 and 48 ns in the test data), while the events of a telescope are hundreds of µs
    apart at the usual trigger rates: the default ``max_time_diff`` (10 µs) does not give random
    coincidences.

    Parameters
    ----------
    config : traitlets.config.Config, optional
        Configuration of the combiner.
    parent : ctapipe.core.Component or ctapipe.core.Tool, optional
        Parent of the combiner, whose configuration is used.
    **kwargs
        Values of the configuration options, e.g. ``max_time_diff=5 * u.us``.
    """

    max_time_diff = AstroQuantity(
        default_value=10 * u.us,
        physical_type=u.physical.time,
        help="Maximum difference of the trigger times of coincident events (strict)",
    ).tag(config=True)

    def match_pair(self, reference, other):
        """
        Mutual closest events of the two telescopes, within ``max_time_diff``.

        Parameters
        ----------
        reference : astropy.table.Table
            Events of the reference telescope, sorted by trigger time (column ``time_ns``).
        other : astropy.table.Table
            Events of the other telescope, sorted by trigger time (column ``time_ns``).

        Returns
        -------
        reference_index, other_index : numpy.ndarray
            Rows of the coincident events in each table.
        """
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
    Coincident events: trigger times within a window, after the search of the offset between the clocks.

    For data without White Rabbit timestamps (taken before 25 September 2023), whose clocks have an
    offset. An event of the other telescope is coincident with an event of the reference telescope
    if its trigger time is within ``window_half_width`` (strictly) of the trigger time of the event
    of the reference telescope plus the offset of the clocks (time of the other telescope -
    time of the reference telescope).

    The offset is ``time_offset`` if given. Otherwise, it is the offset giving the most coincident
    events among the offsets from ``offset_start`` to ``offset_stop`` (excluded) in steps of
    ``offset_step`` (the first one if several give as many), kept in `last_time_offset`.

    An event with several events of the other telescope in its window is ambiguous and not
    matched, in both directions: an event of the reference telescope is matched only if its window
    has a single event of the other telescope, whose window (time - offset) has a single event of
    the reference telescope.

    With the default window (1 ms), the method gives random coincidences at the usual trigger
    rates: it is meant for data without precise timestamps, otherwise use
    `WhiteRabbitClosestCombiner` (or a fixed ``time_offset`` with a narrow window).

    Parameters
    ----------
    config : traitlets.config.Config, optional
        Configuration of the combiner.
    parent : ctapipe.core.Component or ctapipe.core.Tool, optional
        Parent of the combiner, whose configuration is used.
    **kwargs
        Values of the configuration options, e.g. ``time_offset=0 * u.s``.

    Attributes
    ----------
    last_time_offset : int or None
        Offset of the clocks (ns) used by the last call of `match_pair`, None before.
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
        self.last_time_offset = None

    @staticmethod
    def _window(sorted_times, times, half_width):
        """first index and number of ``sorted_times`` strictly within ``half_width`` of each of ``times``"""
        start = np.searchsorted(sorted_times, times - half_width, side="right")
        stop = np.searchsorted(sorted_times, times + half_width, side="left")
        return start, stop - start

    def n_coincidences(self, time_1, time_2, offset_ns):
        """
        Number of coincidences for an offset of the clocks.

        Parameters
        ----------
        time_1 : numpy.ndarray
            Trigger times of the reference telescope (integer ns).
        time_2 : numpy.ndarray
            Trigger times of the other telescope (integer ns), sorted.
        offset_ns : int
            Offset of the clocks (ns): time of the other telescope - time of the reference telescope.

        Returns
        -------
        int
            Number of events of ``time_2`` within ``window_half_width`` (strictly) of
            ``time_1 + offset_ns``, summed over the events of ``time_1``.
        """
        return int(self._window(time_2, time_1 + offset_ns, _ns(self.window_half_width))[1].sum())

    def find_time_offset(self, time_1, time_2):
        """
        Offset of the clocks giving the most coincidences.

        Parameters
        ----------
        time_1 : numpy.ndarray
            Trigger times of the reference telescope (integer ns).
        time_2 : numpy.ndarray
            Trigger times of the other telescope (integer ns), sorted.

        Returns
        -------
        int
            Offset (ns), from ``offset_start`` to ``offset_stop`` (excluded) in steps of
            ``offset_step``, with the most coincidences (`n_coincidences`): the first one if
            several give as many.

        Raises
        ------
        ValueError
            If there is no offset between ``offset_start`` and ``offset_stop``.
        """
        offsets = np.arange(_ns(self.offset_start), _ns(self.offset_stop), _ns(self.offset_step), dtype=np.int64)
        if len(offsets) == 0:
            raise ValueError("SlidingWindowCombiner: no offset between offset_start and offset_stop")
        counts = [self.n_coincidences(time_1, time_2, offset) for offset in offsets]
        best = int(np.argmax(counts))
        self.log.info("Time offset found: %d ns (%d coincidences)", offsets[best], counts[best])
        return int(offsets[best])

    def match_pair(self, reference, other):
        """
        Single events of the other telescope within the window of the events of the reference telescope.

        Parameters
        ----------
        reference : astropy.table.Table
            Events of the reference telescope, sorted by trigger time (column ``time_ns``).
        other : astropy.table.Table
            Events of the other telescope, sorted by trigger time (column ``time_ns``).

        Returns
        -------
        reference_index, other_index : numpy.ndarray
            Rows of the coincident events in each table.

        Raises
        ------
        ValueError
            If the offset is searched and there is no offset between ``offset_start`` and ``offset_stop``.
        """
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
    Coincident events: same array event id given by SWAT.

    For data taken after 30 January 2024: SWAT gives the coincident events of the telescopes the
    same array event id (``arrayEvtNum`` of the raw data), which is the ``event_id`` of the
    events in the DL1 files. Only the event ids are used: the events without the same id are not
    matched, whatever their trigger times, and nothing else is tried if the files have no SWAT
    event ids. An id found several times in the events of a telescope is not matched.

    Parameters
    ----------
    config : traitlets.config.Config, optional
        Configuration of the combiner.
    parent : ctapipe.core.Component or ctapipe.core.Tool, optional
        Parent of the combiner, whose configuration is used.
    **kwargs
        Values of the configuration options, e.g. ``reference_tel_id=22``.
    """

    def match_pair(self, reference, other):
        """
        Events of the two telescopes with the same event id.

        Parameters
        ----------
        reference : astropy.table.Table
            Events of the reference telescope (column ``event_id``).
        other : astropy.table.Table
            Events of the other telescope (column ``event_id``).

        Returns
        -------
        reference_index, other_index : numpy.ndarray
            Rows of the coincident events in each table, without the ids found several times in
            a table.
        """
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
