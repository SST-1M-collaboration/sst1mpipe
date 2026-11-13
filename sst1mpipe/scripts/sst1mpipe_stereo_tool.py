"""
Stereo events of the SST-1M telescopes, from the DL1 files of each telescope (sst1mpipe-process):
the coincident events of the telescopes are combined into stereo events.

The telescopes have no array trigger: each telescope writes its own events (one DL1 file per run
and per telescope). The events of each telescope are read from its DL1 files (any number of files
of each telescope, e.g. the runs of a night), the coincident events are found from the trigger
tables of the files and combined into one event with a StereoEventCombiner, chosen by
``--combiner``:

- WhiteRabbitClosestCombiner (default): closest White Rabbit trigger times, within 10 us
- SlidingWindowCombiner: trigger times within a window, after the search of the offset between
  the clocks of the telescopes (data without White Rabbit timestamps)
- SWATEventIdCombiner: same array event id given by SWAT

A stereo event is an event of the reference telescope (``--reference-tel``, 21 by default) with
the coincident events of the other telescopes, and has the obs_id and event_id of the event of
the reference telescope. Each event is in at most one stereo event. Only the events triggered by the
cameras (SUBARRAY) are combined, the stereo events are SUBARRAY events at the earliest trigger
time of their telescopes. The other events (mono events) are written as they are with
``--keep-mono``, otherwise not. The events are written in time order, with the ctapipe DataWriter
(which data levels are written: ``DataWriter.write_dl1_images``, ``write_dl1_parameters``, ...).

The stereo reconstruction (e.g. the geometry of the showers with the ctapipe ShowerProcessor) is
not done: the stereo DL2 of the events is empty. The output file has the subarray of all the
telescopes, the observation and scheduling blocks of the input files, and their monitoring data
(pointing, statistics of the pedestal events) with ``--StereoTool.copy_monitoring`` (default). The
summary of the processing of the input files (/dl1/info) is not copied.

Usage::

    sst1mpipe-stereo --input SST1M1_20260120_1179.dl1.h5 --input SST1M2_20260120_1102.dl1.h5 \\
        --output SST1M_20260120_1179.stereo.dl1.h5 [--combiner SWATEventIdCombiner] [--keep-mono]
"""
from collections import Counter

import numpy as np
from astropy.table import vstack
from ctapipe.core import Tool, ToolConfigurationError
from ctapipe.core.traits import Bool, ComponentName, List, Path, classes_with_traits, flag
from ctapipe.instrument import SubarrayDescription
from ctapipe.io import DataWriter, HDF5EventSource, read_table
from tqdm import tqdm

from sst1mpipe.io import append_monitoring_tables
from sst1mpipe.stereo import StereoEventCombiner, to_sst1m_event
from sst1mpipe.utils import get_trigger_time_ns, merge_subarrays

__all__ = ["StereoTool", "CombinedRuns", "TelescopeEventReader", "read_telescope_triggers"]

TELESCOPE_TRIGGER_TABLE = "/dl1/event/telescope/trigger"


def read_telescope_triggers(files):
    """
    Trigger table of the events of a telescope (``/dl1/event/telescope/trigger``) in its DL1 files,
    with the column ``file_index`` (index of the file in ``files``) and ``time_ns`` (trigger time, ns)
    """
    triggers = []
    for file_index, path in enumerate(files):
        table = read_table(path, TELESCOPE_TRIGGER_TABLE)
        table["file_index"] = np.full(len(table), file_index)
        triggers.append(table)
    table = vstack(triggers, metadata_conflicts="silent")
    table["time_ns"] = get_trigger_time_ns({"trigger_time": table["time"]})
    return table


class CombinedRuns:
    """
    Run information of the input files of all the telescopes, given to the DataWriter in place of
    an event source: subarray of the telescopes, observation and scheduling blocks
    """

    def __init__(self, sources, subarray):
        self.subarray = subarray
        self.is_simulation = False
        self.observation_blocks, self.scheduling_blocks = {}, {}
        for source in sources:
            for obs_id, block in source.observation_blocks.items():
                self.observation_blocks.setdefault(obs_id, block)
            for sb_id, block in source.scheduling_blocks.items():
                self.scheduling_blocks.setdefault(sb_id, block)
        self.obs_ids = sorted(self.observation_blocks)
        self.simulation_config = {}
        self.atmosphere_density_profile = None


class TelescopeEventReader:
    """
    Events of a telescope read from its DL1 files, by (file index, obs_id, event_id): the files are
    read in order, the events read before they are asked for are kept until they are.
    """

    def __init__(self, sources, wanted):
        self._sources = sources
        self._wanted = set(wanted)
        self._cache = {}
        self._events = self._read()

    def _read(self):
        for file_index, source in enumerate(self._sources):
            for event in source:
                yield (file_index, event.index.obs_id, event.index.event_id), event

    def get(self, key):
        """Event ``key`` (file index, obs_id, event_id)"""
        while key not in self._cache:
            try:
                read_key, event = next(self._events)
            except StopIteration:
                raise KeyError(f"Event {key} not found in the files of the telescope") from None
            if read_key in self._wanted:
                self._cache[read_key] = event
        return self._cache.pop(key)


class StereoTool(Tool):
    """Combination of the coincident events of the SST-1M telescopes into stereo events"""

    name = "sst1mpipe-stereo"
    description = __doc__
    examples = """
    sst1mpipe-stereo -i SST1M1_20260120_1179.dl1.h5 -i SST1M2_20260120_1102.dl1.h5 -o stereo.dl1.h5
    """

    input_files = List(
        Path(exists=True, directory_ok=False, file_ok=True),
        default_value=[],
        help="DL1 files of the telescopes (one telescope per file), any number of files per telescope",
    ).tag(config=True)
    combiner_type = ComponentName(
        StereoEventCombiner,
        default_value="WhiteRabbitClosestCombiner",
        help="Method of the matching of the coincident events of the telescopes",
    ).tag(config=True)
    keep_mono = Bool(
        False, help="Write the events which are not combined (mono events), as they are",
    ).tag(config=True)
    copy_monitoring = Bool(
        True, help="Copy the monitoring data of the input files (pointing, statistics of the pedestal events)",
    ).tag(config=True)
    progress_bar = Bool(False, help="Show a progress bar").tag(config=True)

    aliases = {
        ("i", "input"): "StereoTool.input_files",
        ("o", "output"): "DataWriter.output_path",
        "combiner": "StereoTool.combiner_type",
        "reference-tel": "StereoEventCombiner.reference_tel_id",
    }
    flags = {
        **flag("keep-mono", "StereoTool.keep_mono", "write the mono events", "do not write the mono events"),
        **flag("progress", "StereoTool.progress_bar", "show a progress bar", "do not show a progress bar"),
        **flag("overwrite", "DataWriter.overwrite", "overwrite the output file", "do not overwrite the output file"),
    }
    classes = [DataWriter] + classes_with_traits(StereoEventCombiner)

    def setup(self):
        if len(self.input_files) == 0:
            raise ToolConfigurationError("No input files (--input)")
        self.combiner = StereoEventCombiner.from_name(self.combiner_type, parent=self)

        # input files by telescope
        self.files = {}
        subarrays = {}
        for path in self.input_files:
            subarray = SubarrayDescription.from_hdf(path)
            if len(subarray.tel_ids) != 1:
                raise ToolConfigurationError(
                    f"{path}: input files of one telescope only, got telescopes {subarray.tel_ids.tolist()}"
                )
            tel_id = int(subarray.tel_ids[0])
            self.files.setdefault(tel_id, []).append(path)
            subarrays.setdefault(tel_id, subarray)
        if len(self.files) < 2:
            raise ToolConfigurationError(f"Input files of at least two telescopes needed, got telescopes {sorted(self.files)}")
        if self.combiner.reference_tel_id not in self.files:
            raise ToolConfigurationError(
                f"No input files of the reference telescope {self.combiner.reference_tel_id}"
                f" (--reference-tel), telescopes: {sorted(self.files)}"
            )

        # trigger tables of the telescopes; the files of each telescope in time order
        self.triggers = {}
        for tel_id, files in self.files.items():
            triggers = read_telescope_triggers(files)
            first_times = [
                np.min(triggers["time_ns"][triggers["file_index"] == i], initial=np.iinfo(np.int64).max)
                for i in range(len(files))
            ]
            order = np.argsort(first_times, kind="stable")
            self.files[tel_id] = [files[i] for i in order]
            triggers["file_index"] = np.argsort(order)[triggers["file_index"]]
            self.triggers[tel_id] = triggers[triggers["tel_id"] == tel_id]

        # the HDF5EventSource interpolates the pointing of the telescopes at the time of each event
        # (monitoring data of the file), not written in the output file: an event out of the time
        # range of the pointing table has a NaN pointing instead of stopping the tool
        if "bounds_error" not in self.config.get("PointingInterpolator", {}):
            self.config.PointingInterpolator.bounds_error = False
        self.sources = {
            tel_id: [self.enter_context(HDF5EventSource(path, parent=self)) for path in files]
            for tel_id, files in self.files.items()
        }
        all_sources = [source for sources in self.sources.values() for source in sources]
        if any(source.is_simulation for source in all_sources):
            raise ToolConfigurationError(
                "Simulated input files: the stereo events of the simulations are combined by sim_telarray"
            )
        subarray = merge_subarrays([subarrays[tel_id] for tel_id in sorted(subarrays)])
        self.writer = DataWriter(event_source=CombinedRuns(all_sources, subarray), parent=self)
        # closed in finish(), or when the tool exits if it fails before
        self._writer_closed = False
        self._exit_stack.callback(self._close_writer)
        self.counts = Counter()
        self.used_files = set()

    def _close_writer(self):
        if not self._writer_closed:
            self._writer_closed = True
            self.writer.finish()

    def _output_events(self):
        """
        The events written: by telescope, the row of its event in its trigger table, in time order;
        the stereo events, then the mono events with keep_mono
        """
        groups = self.combiner.match(self.triggers)
        n_stereo = len(groups[self.combiner.reference_tel_id])
        events = [
            {tel_id: rows[i] for tel_id, rows in groups.items() if rows[i] >= 0}
            for i in range(n_stereo)
        ]
        self.counts["stereo"] = n_stereo
        for tel_id, triggers in self.triggers.items():
            used = np.zeros(len(triggers), dtype=bool)
            used[groups[tel_id][groups[tel_id] >= 0]] = True
            self.counts[f"mono_tel_{tel_id:03d}"] = np.count_nonzero(~used)
            if self.keep_mono:
                events.extend({tel_id: row} for row in np.flatnonzero(~used))

        def time(event):
            return min(self.triggers[tel_id]["time_ns"][row] for tel_id, row in event.items())

        events.sort(key=time)
        self._check_unique_ids(events)
        return events

    def _event_id(self, event):
        """(obs_id, event_id) of an output event: those of the reference telescope for a stereo event"""
        tel_id = self.combiner.reference_tel_id if len(event) > 1 else next(iter(event))
        row = self.triggers[tel_id][event[tel_id]]
        return int(row["obs_id"]), int(row["event_id"])

    def _check_unique_ids(self, events):
        ids = Counter(self._event_id(event) for event in events)
        duplicated = [event_id for event_id, n in ids.items() if n > 1]
        if duplicated:
            raise ToolConfigurationError(
                f"{len(duplicated)} output events with the same obs_id and event_id, e.g. {duplicated[0]}:"
                " the mono events of the other telescopes have the ids of events of the reference telescope"
            )

    def _key(self, tel_id, row):
        trigger = self.triggers[tel_id][row]
        return int(trigger["file_index"]), int(trigger["obs_id"]), int(trigger["event_id"])

    def start(self):
        events = self._output_events()
        readers = {}
        for tel_id in self.files:
            wanted = [self._key(tel_id, event[tel_id]) for event in events if tel_id in event]
            readers[tel_id] = TelescopeEventReader(self.sources[tel_id], wanted)
        for count, rows in enumerate(tqdm(events, desc="Stereo events", disable=not self.progress_bar)):
            telescope_events = {}
            for tel_id, row in rows.items():
                key = self._key(tel_id, row)
                telescope_events[tel_id] = readers[tel_id].get(key)
                self.used_files.add(self.files[tel_id][key[0]])
            if len(telescope_events) > 1:
                event = self.combiner.combine(telescope_events)
            else:
                event = to_sst1m_event(next(iter(telescope_events.values())))
            event.count = count
            self.writer(event)

    def finish(self):
        self._close_writer()
        if self.copy_monitoring:
            for path in [path for files in self.files.values() for path in files if path in self.used_files]:
                append_monitoring_tables(path, self.writer.output_path)
        self.log.info("%d stereo events", self.counts["stereo"])
        for name, n in sorted(self.counts.items()):
            if name.startswith("mono"):
                self.log.info(
                    "%d mono events of telescope %s %s", n, name.split("_")[-1],
                    "written" if self.keep_mono else "not written (--keep-mono to write them)",
                )


def main():
    tool = StereoTool()
    tool.run()


if __name__ == "__main__":
    main()
