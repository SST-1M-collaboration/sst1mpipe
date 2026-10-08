from ctapipe.io import EventSource
import threading
from itertools import islice

import numpy as np
import tables
import zmq
import pytest
from ctapipe.core import run_tool
from ctapipe.io import read_table
from traitlets import traitlets


from sst1mpipe.io.zmq_event_source import ZMQEventSource, run_header_obs_id
from sst1mpipe.constants import N_PIXELS, N_CHANNELS, SUBARRAY_DESCRIPTION
from ctapipe.containers import EventType, SchedulingBlockType
from ctapipe.io.datalevels import DataLevel
from protozfits import (
    CoreMessages_pb2,
    DL0v1_Telescope_pb2,
    File,
    ProtoDataModel_pb2,
    R1v1_pb2,
    make_namedtuple,
    numpy_to_any_array,
)

from sst1mpipe.io.containers import StreamCameraConfigContainer
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.resources import RTA_CONFIG_FILE, TEST_DATA_DIR
from sst1mpipe.scripts.sst1mpipe_process_tool import ProcessorTool

def create_fake_dl0_event_message(event_id: int, tel_id: int ):

    event = DL0v1_Telescope_pb2.Event()
    event.event_id = event_id
    event.tel_id = tel_id
    event.event_type = 32
    event.event_time_s = event_id
    event.event_time_qns = 4
    event.num_channels = N_CHANNELS
    event.num_samples = 50
    event.num_pixels_survived = N_PIXELS
    event.waveform.CopyFrom(numpy_to_any_array(np.ones((N_PIXELS, event.num_samples))))
    event.pedestal_intensity.CopyFrom(numpy_to_any_array(np.zeros(N_PIXELS)))

    message = CoreMessages_pb2.CTAMessage()
    message.payload_type.append(
        CoreMessages_pb2.DL0_TELESCOPE_EVENT
    )
    message.payload_data.append(
        event.SerializeToString()
    )
    message.source_name = "pytest"

    return message.SerializeToString()

@pytest.mark.parametrize("n_events", [1000])
def test_zmq_event_source(n_events):

    endpoint = "inproc://test"
    source = ZMQEventSource(endpoint, max_events=n_events)
    tel_id = 1

    producer = source.socket.context.socket(zmq.PUSH)
    producer.bind(endpoint)

    for i in range(n_events):
        producer.send(create_fake_dl0_event_message(i, tel_id))

    k = 0
    for event in source:


        assert event.count == k
        assert event.dl0.tel[tel_id].waveform.sum() == N_CHANNELS * N_PIXELS * 50
        assert event.index.event_id == k
        k += 1

    assert event.count == n_events - 1

@pytest.mark.parametrize(["endpoint", "valid"],
                         [  ["inproc://test", True],
                            ["tcp://192.168.1.1:1986", True],
                            ["tcp://[2a7d:91c4:8f21:3b7a:5e12:aa90:1c44:72ef]:8000", True],
                            [ "not_a_valid_endpoint", False],
                            [ "/some/folder/on/linux", False],
                          ])
def test_zmq_address(endpoint, valid):

    assert ZMQEventSource.is_compatible(endpoint) == valid

def test_subarray(tmp_path):

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    file = data_dir / "dummy-array.h5"

    subarray = SUBARRAY_DESCRIPTION
    subarray.to_hdf(file)
    source = ZMQEventSource(input_url="inproc://test", subarray_file=file)

    assert subarray.name == source.subarray.name

@pytest.mark.xfail(raises=traitlets.TraitError) # Unfortunately ctapipe does not allow urls with tcp:// but only Path
@pytest.mark.parametrize("url", ["tcp://192.168.1.1:1986", "inproc://test"])
def test_zmq_event_source_from_event_source(url):

    EventSource(input_url=url)


# ---------------------------------------------------------------------------
# stream with the data stream, the camera configuration and the end of stream
# ---------------------------------------------------------------------------

SB_ID, OBS_ID, LOCAL_RUN_ID = 1000, 2000, 33
WAVEFORM_SCALE, WAVEFORM_OFFSET = 4.0, 10.0
N_SAMPLES = 50


def cta_message(*payloads):
    """CTAMessage with the (type, protobuf message) payloads"""
    message = CoreMessages_pb2.CTAMessage()
    for msg_type, payload in payloads:
        message.payload_type.append(msg_type)
        message.payload_data.append(b"" if payload is None else payload.SerializeToString())
    message.source_name = "pytest"
    return message.SerializeToString()


def r1_data_stream(tel_id):
    stream = R1v1_pb2.TelescopeDataStream()
    stream.tel_id, stream.sb_id, stream.obs_id = tel_id, SB_ID, OBS_ID
    stream.waveform_scale, stream.waveform_offset = WAVEFORM_SCALE, WAVEFORM_OFFSET
    return CoreMessages_pb2.TELESCOPE_DATA_STREAM, stream


def r1_camera_config(tel_id):
    config = R1v1_pb2.CameraConfiguration()
    config.tel_id, config.local_run_id, config.config_time_s = tel_id, LOCAL_RUN_ID, 1769015262.5
    config.camera_config_id, config.num_modules, config.num_pixels = 7, 108, N_PIXELS
    config.num_channels, config.num_samples_nominal = N_CHANNELS, N_SAMPLES
    config.data_model_version = "1.0"
    config.pixel_id_map.CopyFrom(numpy_to_any_array(np.arange(N_PIXELS, dtype=np.uint16)))
    config.module_id_map.CopyFrom(numpy_to_any_array(np.arange(108, dtype=np.uint16)))
    return CoreMessages_pb2.CAMERA_CONFIG, config


def r1_event(event_id, tel_id, stored_value, event_type=32):
    event = R1v1_pb2.Event()
    event.event_id, event.tel_id, event.local_run_id, event.event_type = event_id, tel_id, LOCAL_RUN_ID, event_type
    event.event_time_s, event.event_time_qns = 1769015262 + event_id, 4 * 383583536
    event.num_channels, event.num_pixels, event.num_samples = N_CHANNELS, N_PIXELS, N_SAMPLES
    event.waveform.CopyFrom(numpy_to_any_array(np.full((N_PIXELS, N_SAMPLES), stored_value, dtype=np.uint16)))
    event.pedestal_intensity.CopyFrom(numpy_to_any_array(np.zeros(N_PIXELS, dtype=np.float32)))
    event.pixel_status.CopyFrom(numpy_to_any_array(np.ones(N_PIXELS, dtype=np.uint8)))
    return CoreMessages_pb2.R1_EVENT, event


END_OF_STREAM = (CoreMessages_pb2.END_OF_STREAM, None)


def send(source, messages, endpoint):
    producer = source.socket.context.socket(zmq.PUSH)
    producer.bind(endpoint)
    for message in messages:
        producer.send(message)
    return producer


def test_r1_stream():
    endpoint = "inproc://test_r1_stream"
    tel_id = 21
    source = ZMQEventSource(endpoint)
    send(source, [
        cta_message(r1_data_stream(tel_id), r1_camera_config(tel_id)),
        # several events in a message
        cta_message(r1_event(1, tel_id, 50), r1_event(2, tel_id, 90, event_type=2)),
        cta_message(r1_event(3, tel_id, 130, event_type=200)),
        cta_message(END_OF_STREAM),
    ], endpoint)

    # the reading stops at the end of the stream
    events = list(source)
    assert [event.count for event in events] == [0, 1, 2]
    assert [event.index.event_id for event in events] == [1, 2, 3]
    assert {event.index.obs_id for event in events} == {OBS_ID}

    # waveform decoded with the scale and offset of the data stream
    for event, stored in zip(events, [50, 90, 130], strict=True):
        r1 = event.r1.tel[tel_id]
        assert r1.waveform.shape == (N_CHANNELS, N_PIXELS, N_SAMPLES)
        assert r1.waveform.dtype == np.float32
        np.testing.assert_allclose(r1.waveform, stored / WAVEFORM_SCALE - WAVEFORM_OFFSET)
        # a container for each event
        assert list(event.r1.tel) == [tel_id]

    # trigger
    assert [event.trigger.event_type for event in events] == [EventType.SUBARRAY, EventType.SKY_PEDESTAL, EventType.UNKNOWN]
    for event in events:
        assert event.trigger.tels_with_trigger == [tel_id]
        assert event.trigger.time == event.r1.tel[tel_id].event_time
        assert event.trigger.tel[tel_id].time == event.trigger.time
    assert events[0].trigger.time.unix_tai == pytest.approx(1769015263.383583536, abs=1e-6)

    # blocks of the data stream
    assert list(source.scheduling_blocks) == [SB_ID]
    assert source.scheduling_blocks[SB_ID].sb_type == SchedulingBlockType.OBSERVATION
    assert source.scheduling_blocks[SB_ID].producer_id == "SST1M-21"
    assert list(source.observation_blocks) == [OBS_ID]
    assert source.observation_blocks[OBS_ID].sb_id == SB_ID
    assert source.data_streams[tel_id].waveform_scale == WAVEFORM_SCALE

    # camera configuration
    config = source.camera_config[tel_id]
    assert isinstance(config, StreamCameraConfigContainer)
    assert config.data_level == "R1"
    assert (config.local_run_id, config.num_pixels, config.num_samples_nominal) == (LOCAL_RUN_ID, N_PIXELS, N_SAMPLES)
    assert config.config_time_s == 1769015262.5
    np.testing.assert_array_equal(config.pixel_id_map, np.arange(N_PIXELS))
    source.close()


def test_r1_event_before_data_stream():
    """without data stream: the stored waveform and the local run id of the camera as obs_id"""
    endpoint = "inproc://test_r1_event_before_data_stream"
    source = ZMQEventSource(endpoint)
    send(source, [cta_message(r1_event(1, 22, 50)), cta_message(END_OF_STREAM)], endpoint)

    (event,) = list(source)
    assert event.index.obs_id == LOCAL_RUN_ID
    np.testing.assert_allclose(event.r1.tel[22].waveform, 50)
    assert source.scheduling_blocks == {}
    source.close()


def test_dl0_stream():
    endpoint = "inproc://test_dl0_stream"
    tel_id = 22
    stream = DL0v1_Telescope_pb2.DataStream()
    stream.tel_id, stream.sb_id, stream.obs_id = tel_id, SB_ID, OBS_ID
    stream.waveform_scale, stream.waveform_offset = WAVEFORM_SCALE, WAVEFORM_OFFSET
    config = DL0v1_Telescope_pb2.CameraConfiguration()
    config.tel_id, config.local_run_id, config.num_pixels, config.sampling_frequency = tel_id, LOCAL_RUN_ID, N_PIXELS, 250

    source = ZMQEventSource(endpoint, scheduling_block_type="CALIBRATION")
    message = CoreMessages_pb2.CTAMessage()
    message.ParseFromString(create_fake_dl0_event_message(5, tel_id))
    send(source, [
        cta_message((CoreMessages_pb2.DL0_TELESCOPE_DATA_STREAM, stream), (CoreMessages_pb2.DL0_TELESCOPE_CAMERA_CONFIG, config)),
        message.SerializeToString(),
        cta_message(END_OF_STREAM),
    ], endpoint)

    (event,) = list(source)
    assert event.index.obs_id == OBS_ID
    assert event.trigger.event_type == EventType.SUBARRAY
    # stored waveform of 1, pedestal of 0
    np.testing.assert_allclose(event.dl0.tel[tel_id].waveform, 1 / WAVEFORM_SCALE - WAVEFORM_OFFSET)
    assert source.scheduling_blocks[SB_ID].sb_type == SchedulingBlockType.CALIBRATION
    assert source.camera_config[tel_id].data_level == "DL0"
    assert source.camera_config[tel_id].sampling_frequency == 250
    source.close()


def test_process_r1_stream(tmp_path):
    """sst1mpipe-process reading a R1 stream up to the DL1 parameters"""
    tel_id, n_events = 21, 20
    context = zmq.Context()
    producer = context.socket(zmq.PUSH)
    port = producer.bind_to_random_port("tcp://127.0.0.1")
    rng = np.random.default_rng(0)
    events = []
    for i in range(n_events):
        message_type, event = r1_event(i, tel_id, 0)
        # stored values of a waveform of 0 +- 1 (decoded), with a signal in a few pixels
        waveform = rng.normal(0, 1, (N_PIXELS, N_SAMPLES))
        waveform[:20, 20:25] += 30
        stored = np.round((waveform + WAVEFORM_OFFSET) * WAVEFORM_SCALE).astype(np.uint16)
        event.waveform.CopyFrom(numpy_to_any_array(stored))
        events.append(cta_message((message_type, event)))
    messages = [cta_message(r1_data_stream(tel_id), r1_camera_config(tel_id)), *events, cta_message(END_OF_STREAM)]
    # the messages are sent once the tool is connected
    sender = threading.Thread(target=lambda: [producer.send(message) for message in messages], daemon=True)
    sender.start()

    output = tmp_path / "events.dl1.h5"
    try:
        run_tool(ProcessorTool(), argv=[
            f"--input=tcp://127.0.0.1:{port}",
            f"--output={output}",
            f"--config={RTA_CONFIG_FILE}",
        ], raises=True)
    finally:
        sender.join(timeout=10)
        producer.close(linger=0)
        context.term()

    trigger = read_table(output, "/dl1/event/subarray/trigger")
    assert len(trigger) == n_events
    assert set(trigger["obs_id"]) == {OBS_ID}
    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
    assert len(parameters) == n_events
    # the blocks are written by the DataWriter before the data stream message is received
    with tables.open_file(output) as h5:
        assert "/configuration/observation/observation_block" not in h5


# ---------------------------------------------------------------------------
# DigiCam camera events (DataModel.CameraEvent), replayed from a zfits file
# ---------------------------------------------------------------------------

ZFITS_FILE = TEST_DATA_DIR / "zfits" / "SST1M1_20260121_0206.fits.fz"
ZFITS_TEL_ID = 21
R0_FIELDS = [
    "waveform", "pedestal", "pixel_flags", "camera_event_number", "event_type",
    "trigger_input_traces", "trigger_output_patch7", "trigger_output_patch19", "trigger_output_muon",
]


def camera_event_messages(path, n_events=None, modify=None):
    """CAMERA_EVENT messages of the events (protobuf) of a zfits file"""
    with File(str(path), pure_protobuf=True) as f:
        events = list(islice(f.Events, n_events))
    messages = []
    for event in events:
        if modify is not None:
            modify(event)
        message = CoreMessages_pb2.CTAMessage()
        message.payload_type.append(CoreMessages_pb2.CAMERA_EVENT)
        message.payload_data.append(event.SerializeToString())
        messages.append(message.SerializeToString())
    return messages


def assert_same_r0_events(stream_events, file_source):
    file_events = iter(file_source)
    for stream_event in stream_events:
        file_event = next(file_events)
        assert stream_event.index.event_id == file_event.index.event_id
        assert stream_event.trigger.time == file_event.trigger.time
        assert stream_event.trigger.event_type == file_event.trigger.event_type
        assert stream_event.trigger.tels_with_trigger == file_event.trigger.tels_with_trigger
        for name in R0_FIELDS:
            np.testing.assert_array_equal(
                getattr(stream_event.r0.tel[ZFITS_TEL_ID], name), getattr(file_event.r0.tel[ZFITS_TEL_ID], name),
            )
    assert next(file_events, None) is None


def test_camera_event_stream():
    """the camera events of the stream are read as the events of the zfits file"""
    endpoint = "inproc://test_camera_event_stream"
    source = ZMQEventSource(endpoint)
    send(source, [
        cta_message(r1_data_stream(ZFITS_TEL_ID)),
        *camera_event_messages(ZFITS_FILE),
        cta_message(END_OF_STREAM),
    ], endpoint)

    events = list(source)
    assert [event.count for event in events] == list(range(130))
    assert {event.index.obs_id for event in events} == {OBS_ID}
    assert DataLevel.R0 in source.datalevels
    # SWAT array event ids, as for the file
    assert_same_r0_events(events, SST1MEventSource(ZFITS_FILE))
    source.close()


def test_camera_event_stream_swapped_modules():
    swapped_modules = [{"tel_id": ZFITS_TEL_ID, "start": "2026-01-21T00:00:00", "stop": "2026-01-22T00:00:00", "modules": [8, 9]}]
    endpoint = "inproc://test_camera_event_stream_swapped_modules"
    source = ZMQEventSource(endpoint, swapped_modules=swapped_modules)
    send(source, [*camera_event_messages(ZFITS_FILE, n_events=5), cta_message(END_OF_STREAM)], endpoint)

    events = list(source)
    assert_same_r0_events(events, SST1MEventSource(ZFITS_FILE, max_events=5, swapped_modules=swapped_modules))
    # the waveforms are not those without swap
    unswapped = next(iter(SST1MEventSource(ZFITS_FILE, max_events=1)))
    assert not np.array_equal(events[0].r0.tel[ZFITS_TEL_ID].waveform, unswapped.r0.tel[ZFITS_TEL_ID].waveform)
    source.close()


def test_camera_event_ids_without_swat():
    """without SWAT array event id in the first event, the camera event number is the event id"""
    def remove_swat_id(event):
        event.arrayEvtNum = 0

    endpoint = "inproc://test_camera_event_ids_without_swat"
    source = ZMQEventSource(endpoint)
    send(source, [*camera_event_messages(ZFITS_FILE, n_events=3, modify=remove_swat_id), cta_message(END_OF_STREAM)], endpoint)

    events = list(source)
    assert [event.index.event_id for event in events] == [
        event.r0.tel[ZFITS_TEL_ID].camera_event_number for event in events
    ]
    assert events[0].index.event_id == 58134443
    # no data stream
    assert {event.index.obs_id for event in events} == {0}
    source.close()


def test_process_camera_event_stream(tmp_path):
    """sst1mpipe-process gives the same DL1 parameters for the stream and the zfits file"""
    context = zmq.Context()
    producer = context.socket(zmq.PUSH)
    port = producer.bind_to_random_port("tcp://127.0.0.1")
    messages = [*camera_event_messages(ZFITS_FILE), cta_message(END_OF_STREAM)]
    sender = threading.Thread(target=lambda: [producer.send(message) for message in messages], daemon=True)
    sender.start()

    stream_output = tmp_path / "stream.dl1.h5"
    try:
        tool = ProcessorTool()
        run_tool(tool, argv=[
            f"--input=tcp://127.0.0.1:{port}",
            f"--output={stream_output}",
            f"--config={RTA_CONFIG_FILE}",
        ], raises=True)
    finally:
        sender.join(timeout=10)
        producer.close(linger=0)
        context.term()
    # the R0 -> R1 calibration and the pedestal monitoring of the SST-1M raw data are used
    assert tool.r0_r1_calibrator is not None
    assert tool.r0_pedestal_monitor.n_buffered(ZFITS_TEL_ID) > 0

    file_output = tmp_path / "file.dl1.h5"
    run_tool(ProcessorTool(), argv=[
        f"--input={ZFITS_FILE}",
        f"--output={file_output}",
        f"--config={RTA_CONFIG_FILE}",
        "--ProcessorTool.wobble_in_output_name=False",
        # run of the transition between two wobbles
        "--ProcessorTool.allowed_sb_types=UNKNOWN",
    ], raises=True)

    table = f"/dl1/event/telescope/parameters/tel_{ZFITS_TEL_ID:03d}"
    stream_parameters, file_parameters = read_table(stream_output, table), read_table(file_output, table)
    assert len(stream_parameters) == 130
    np.testing.assert_array_equal(stream_parameters["event_id"], file_parameters["event_id"])
    for column in ["camera_frame_hillas_intensity", "camera_frame_hillas_x", "camera_frame_hillas_width"]:
        np.testing.assert_array_equal(stream_parameters[column], file_parameters[column])
    assert np.isfinite(stream_parameters["camera_frame_hillas_intensity"]).sum() > 0


# ---------------------------------------------------------------------------
# DigiCam run header (DataModel.CameraRunHeader)
# ---------------------------------------------------------------------------

def camera_run_header(tel_id, run_number, date_mjd=0, run_id=0):
    header = ProtoDataModel_pb2.CameraRunHeader()
    header.telescopeID, header.runNumber, header.dateMJD, header.run_id = tel_id, run_number, date_mjd, run_id
    return CoreMessages_pb2.CAMERA_RUN_HEADER, header


@pytest.mark.parametrize("run_number, date_mjd, run_id, expected", [
    (206, 61061, 0, 202601210206),  # MJD 61061: 2026-01-21
    (206, 0, 0, 206),
    (206, 61061, 123456, 123456),
])
def test_run_header_obs_id(run_number, date_mjd, run_id, expected):
    _, header = camera_run_header(21, run_number, date_mjd, run_id)
    assert run_header_obs_id(make_namedtuple(header)) == expected


def test_camera_event_stream_run_header():
    """the obs_id of the camera events is given by the run header of the telescope"""
    endpoint = "inproc://test_camera_event_stream_run_header"
    source = ZMQEventSource(endpoint)
    send(source, [
        cta_message(camera_run_header(ZFITS_TEL_ID, 206, date_mjd=61061)),
        *camera_event_messages(ZFITS_FILE, n_events=3),
        cta_message(END_OF_STREAM),
    ], endpoint)

    events = list(source)
    # as for the zfits file SST1M1_20260121_0206
    assert {event.index.obs_id for event in events} == {202601210206}
    assert source.run_headers[ZFITS_TEL_ID].runNumber == 206
    assert list(source.observation_blocks) == [202601210206]
    assert source.scheduling_blocks[202601210206].sb_type == SchedulingBlockType.OBSERVATION
    assert source.scheduling_blocks[202601210206].producer_id == "SST1M-21"
    source.close()


def test_data_stream_obs_id_before_run_header():
    """the obs_id of the data stream of the telescope is used rather than the one of the run header"""
    endpoint = "inproc://test_data_stream_obs_id_before_run_header"
    source = ZMQEventSource(endpoint)
    send(source, [
        cta_message(camera_run_header(ZFITS_TEL_ID, 206, date_mjd=61061), r1_data_stream(ZFITS_TEL_ID)),
        *camera_event_messages(ZFITS_FILE, n_events=2),
        cta_message(END_OF_STREAM),
    ], endpoint)

    events = list(source)
    assert {event.index.obs_id for event in events} == {OBS_ID}
    assert set(source.observation_blocks) == {202601210206, OBS_ID}
    source.close()
