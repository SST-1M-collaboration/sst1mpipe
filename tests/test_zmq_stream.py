from ctapipe.io import EventSource
from pyvo.io.uws import endpoint
from tqdm import tqdm as tqdm
import numpy as np
import zmq
import pytest
from traitlets import traitlets

from sst1mpipe.io.zmq_event_source import ZMQEventSource
from sst1mpipe.constants import N_PIXELS, N_CHANNELS
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

def test_zmq_event_source():

    endpoint = "inproc://test"
    n_events = int(1E3)
    source = ZMQEventSource(endpoint, max_events=n_events)
    tel_id = 1

    producer = source.socket.context.socket(zmq.PUSH)
    producer.bind(endpoint)

    for i in range(n_events):
        producer.send(create_fake_dl0_event_message(i, tel_id))

    k = 0
    for event in tqdm(source):


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


@pytest.mark.xfail(raises=traitlets.TraitError) # Unfortunately ctapipe does not allow urls but only Path
def test_zmq_event_source_from_event_source():

    EventSource(input_url="tcp://127.0.0.1:5555")
