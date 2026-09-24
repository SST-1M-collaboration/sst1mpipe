from ctapipe.instrument import SubarrayDescription
from ctapipe.io import EventSource
import numpy as np
import zmq
import pytest
from traitlets import traitlets


from sst1mpipe.io.zmq_event_source import ZMQEventSource
from sst1mpipe.constants import N_PIXELS, N_CHANNELS, SUBARRAY_DESCRIPTION
from protozfits import DL0v1_Telescope_pb2, CoreMessages_pb2, numpy_to_any_array

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
