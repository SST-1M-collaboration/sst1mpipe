"""
sst1mpipe-stereo on the observation runs of Mrk 421 (the run of tel 22 is taken during the run of
tel 21): stereo events of their DL1 files (sst1mpipe-process) with the 2 coincident events
COINCIDENT_EVENTS
"""
import shutil

import numpy as np
import pytest
import tables
from ctapipe.containers import EventType
from ctapipe.core import ToolConfigurationError, run_tool
from ctapipe.io import HDF5EventSource, read_table
from traitlets.config import Config

from conftest import COINCIDENT_EVENTS, EVENTS_TEL_1, OBSERVATION_RUNS, run_process
from sst1mpipe.io import MONITORING_GROUPS, attach_monitoring_file
from sst1mpipe.scripts.sst1mpipe_stereo_tool import StereoTool
from sst1mpipe.utils import get_trigger_time_ns

OBS_IDS = {21: 202601201179, 22: 202601201102}
# the monitoring file of tel 21 has the pointing of its first events only (observation_monitoring):
# the pointing of its later events is NaN when the files are read
NAN_POINTING_OUT_OF_RANGE = Config({"PointingInterpolator": {"bounds_error": False}})


@pytest.fixture(scope="module")
def dl1_files(tmp_path_factory, observation_monitoring):
    """
    DL1 files of the observation runs by telescope (the events of tel 21 up to the end of the run
    of tel 22), with the monitoring data of the monitoring files attached
    """
    directory = tmp_path_factory.mktemp("stereo_dl1")
    paths = {}
    for tel_id, run in OBSERVATION_RUNS.items():
        paths[tel_id] = directory / f"tel_{tel_id}.dl1.h5"
        max_events = [f"--max-events={EVENTS_TEL_1.stop}"] if tel_id == 21 else []
        run_process(run["path"], paths[tel_id], *max_events)
        attach_monitoring_file(observation_monitoring[tel_id], paths[tel_id])
    return paths


def run_stereo(inputs, output, *args):
    tool = StereoTool()
    run_tool(tool, argv=[*[f"--input={path}" for path in inputs], f"--output={output}", *args], raises=True)
    return tool


def n_events(path, tel_id):
    return len(read_table(path, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}"))


@pytest.mark.parametrize("combiner_args", [
    ["--combiner=WhiteRabbitClosestCombiner"],
    ["--combiner=SWATEventIdCombiner"],
    # data with White Rabbit timestamps: no clock offset, small window
    ["--combiner=SlidingWindowCombiner", "--SlidingWindowCombiner.time_offset=0 s",
     "--SlidingWindowCombiner.window_half_width=10 us"],
])
def test_stereo_events(dl1_files, tmp_path, combiner_args):
    output = tmp_path / "stereo.dl1.h5"
    tool = run_stereo(dl1_files.values(), output, *combiner_args)
    assert tool.counts["stereo"] == len(COINCIDENT_EVENTS)

    # the coincident events, with the ids of the reference telescope, SUBARRAY events
    trigger = read_table(output, "/dl1/event/subarray/trigger")
    assert trigger["event_id"].tolist() == sorted(COINCIDENT_EVENTS)
    assert np.all(trigger["obs_id"] == OBS_IDS[21])
    assert np.all(trigger["event_type"] == EventType.SUBARRAY.value)
    assert np.all(trigger["tels_with_trigger"])

    # trigger times of the telescopes; stereo event at the earliest one (tel 21 here)
    telescope_trigger = read_table(output, "/dl1/event/telescope/trigger")
    assert len(telescope_trigger) == 2 * len(COINCIDENT_EVENTS)
    times = {
        tel_id: get_trigger_time_ns({"trigger_time": telescope_trigger["time"][telescope_trigger["tel_id"] == tel_id]})
        for tel_id in (21, 22)
    }
    assert (times[22] - times[21]).tolist() == [COINCIDENT_EVENTS[event_id] for event_id in sorted(COINCIDENT_EVENTS)]
    np.testing.assert_array_equal(get_trigger_time_ns({"trigger_time": trigger["time"]}), times[21])
    assert np.all(telescope_trigger["event_type"] == EventType.SUBARRAY.value)
    assert np.all(telescope_trigger["n_trigger_pixels"] > 0)

    # the DL1 parameters of both telescopes, those of their DL1 files
    for tel_id in (21, 22):
        parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
        mono = read_table(dl1_files[tel_id], f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
        mono = mono[np.isin(mono["event_id"], list(COINCIDENT_EVENTS))]
        np.testing.assert_array_equal(parameters["camera_frame_hillas_intensity"], mono["camera_frame_hillas_intensity"])

    with HDF5EventSource(output, config=NAN_POINTING_OUT_OF_RANGE) as source:
        assert source.subarray.tel_ids.tolist() == [21, 22]
        assert sorted(source.obs_ids) == sorted(OBS_IDS.values())
        events = list(source)
    assert [event.index.event_id for event in events] == sorted(COINCIDENT_EVENTS)
    for event in events:
        assert event.trigger.event_type == EventType.SUBARRAY
        assert event.trigger.tels_with_trigger.tolist() == [21, 22]
        assert sorted(event.dl1.tel) == [21, 22]


def test_keep_mono(dl1_files, tmp_path):
    output = tmp_path / "stereo_mono.dl1.h5"
    tool = run_stereo(dl1_files.values(), output, "--keep-mono")
    n_mono = {tel_id: n_events(path, tel_id) - len(COINCIDENT_EVENTS) for tel_id, path in dl1_files.items()}
    assert tool.counts["mono_tel_021"] == n_mono[21] and tool.counts["mono_tel_022"] == n_mono[22]

    trigger = read_table(output, "/dl1/event/subarray/trigger")
    assert len(trigger) == len(COINCIDENT_EVENTS) + n_mono[21] + n_mono[22]
    # in time order, each event once
    assert np.all(np.diff(get_trigger_time_ns({"trigger_time": trigger["time"]})) >= 0)
    assert len(np.unique(trigger[["obs_id", "event_id"]])) == len(trigger)
    n_telescopes = trigger["tels_with_trigger"].sum(axis=1)
    assert np.count_nonzero(n_telescopes == 2) == len(COINCIDENT_EVENTS)
    # the mono events of tel 22 keep the obs_id of their run
    mono_22 = trigger[(n_telescopes == 1) & trigger["tels_with_trigger"][:, 1]]
    assert len(mono_22) == n_mono[22] and np.all(mono_22["obs_id"] == OBS_IDS[22])
    for tel_id in (21, 22):
        assert n_events(output, tel_id) == n_events(dl1_files[tel_id], tel_id)


def test_reference_telescope(dl1_files, tmp_path):
    """with tel 22 as reference: the same stereo events, with the obs_id of the run of tel 22"""
    output = tmp_path / "stereo_22.dl1.h5"
    run_stereo(dl1_files.values(), output, "--reference-tel=22", "--combiner=SWATEventIdCombiner")
    trigger = read_table(output, "/dl1/event/subarray/trigger")
    assert trigger["event_id"].tolist() == sorted(COINCIDENT_EVENTS)
    assert np.all(trigger["obs_id"] == OBS_IDS[22])


def test_monitoring_copied(dl1_files, tmp_path):
    """monitoring data of both telescopes copied, not the summary of the processing (/dl1/info)"""
    output = tmp_path / "stereo.dl1.h5"
    run_stereo(dl1_files.values(), output)

    def monitoring_tables(path):
        with tables.open_file(path) as f:
            return {
                node._v_pathname: node.nrows
                for group in MONITORING_GROUPS if group in f
                for node in f.get_node(group)._f_walknodes("Table")
            }

    copied = monitoring_tables(output)
    for path in dl1_files.values():
        mono = monitoring_tables(path)
        assert len(mono) > 0
        assert all(copied[table] == n_rows for table, n_rows in mono.items())
    assert any("pointing/tel_021" in table for table in copied) and any("pointing/tel_022" in table for table in copied)
    with tables.open_file(output) as f:
        assert "/dl1/info" not in f

    output = tmp_path / "stereo_no_monitoring.dl1.h5"
    run_stereo(dl1_files.values(), output, "--StereoTool.copy_monitoring=False")
    assert monitoring_tables(output) == {}


def test_data_writer_options(dl1_files, tmp_path):
    """the data levels written are those of the DataWriter"""
    output = tmp_path / "stereo.dl1.h5"
    run_stereo(dl1_files.values(), output, "--DataWriter.write_dl1_parameters=False")
    with tables.open_file(output) as f:
        assert "/dl1/event/telescope/parameters" not in f
        assert "/dl1/event/subarray/trigger" in f


def test_configuration_errors(dl1_files, tmp_path):
    output = tmp_path / "stereo.dl1.h5"
    with pytest.raises(ToolConfigurationError, match="two telescopes"):
        run_stereo([dl1_files[21]], output)
    with pytest.raises(ToolConfigurationError, match="reference telescope 23"):
        run_stereo(dl1_files.values(), output, "--reference-tel=23")
    # the events of tel 21 out of the time range of its pointing table stop the tool if asked
    with pytest.raises(ValueError, match="interpolation range"):
        run_stereo(dl1_files.values(), output, "--PointingInterpolator.bounds_error=True", "--overwrite")
    # a file of tel 22 given twice with --keep-mono: its mono events twice
    copy = tmp_path / "tel_22_copy.dl1.h5"
    shutil.copy(dl1_files[22], copy)
    with pytest.raises(ToolConfigurationError, match="same obs_id and event_id"):
        run_stereo([*dl1_files.values(), copy], tmp_path / "duplicated.dl1.h5", "--keep-mono")
