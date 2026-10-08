from ctapipe.io import EventSource
import threading

import numpy as np
import tables
import zmq
import pytest
from ctapipe.core import run_tool
from ctapipe.io import read_table
from traitlets import traitlets


from sst1mpipe.io.zmq_event_source import ZMQEventSource
from sst1mpipe.constants import N_PIXELS, N_CHANNELS, SUBARRAY_DESCRIPTION
from ctapipe.containers import EventType, SchedulingBlockType
from protozfits import DL0v1_Telescope_pb2, CoreMessages_pb2, R1v1_pb2, numpy_to_any_array

from sst1mpipe.io.containers import StreamCameraConfigContainer
from sst1mpipe.resources import RTA_CONFIG_FILE
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
