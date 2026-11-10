"""Stream of the processed events (ZMQEventPublisher -> ZMQEventReceiver)"""

import threading
import time

import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time
from ctapipe.containers import (
    ArrayEventContainer,
    CameraHillasParametersContainer,
    DL1CameraContainer,
    EventType,
    ImageParametersContainer,
    ParticleClassificationContainer,
    ReconstructedEnergyContainer,
    ReconstructedGeometryContainer,
    TelescopeImpactParameterContainer,
)

from sst1mpipe.io.zmq_event_stream import (
    STEREO_TOPIC,
    ZMQEventPublisher,
    ZMQEventReceiver,
    decode,
    encode,
    event_messages,
)


def make_event(event_id=7):
    """A processed event of two telescopes, with DL1, DL2 and a stereo DL2"""
    event = ArrayEventContainer()
    event.index.obs_id, event.index.event_id = 1, event_id
    event.trigger.time = Time(1768964778.5, format="unix")
    event.trigger.event_type = EventType.SUBARRAY
    for tel_id in (21, 22):
        event.trigger.tel[tel_id].time = Time(1768964778.5, format="unix")
        event.r1.tel[tel_id].waveform = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)
        parameters = ImageParametersContainer(hillas=CameraHillasParametersContainer(
            intensity=120.0, x=0.05 * u.m, y=-0.02 * u.m, length=0.03 * u.m, width=0.01 * u.m, psi=30 * u.deg))
        event.dl1.tel[tel_id] = DL1CameraContainer(
            image=np.ones(1296, dtype=np.float32), peak_time=np.zeros(1296, dtype=np.float32),
            image_mask=np.zeros(1296, dtype=bool), is_valid=True, parameters=parameters)
        event.dl2.tel[tel_id].energy["RF"] = ReconstructedEnergyContainer(energy=2 * u.TeV, is_valid=True)
        event.dl2.tel[tel_id].particle_type["Classifier"] = ParticleClassificationContainer(prediction=0.7)
        event.dl2.tel[tel_id].impact["Hillas"] = TelescopeImpactParameterContainer(distance=120 * u.m)
    event.dl2.stereo.geometry["Hillas"] = ReconstructedGeometryContainer(
        alt=1.2 * u.rad, az=20 * u.deg, core_x=10 * u.m, core_y=-5 * u.m, is_valid=True, telescopes=[21, 22])
    return event


def test_event_messages():
    """The fields of the levels have the names of the DataWriter columns, in the units of the tables"""
    messages = event_messages(make_event())
    assert [(topic, tel_id) for topic, tel_id, _ in messages] == [("tel21", 21), ("tel22", 22), (STEREO_TOPIC, 0)]
    levels = messages[0][2]
    assert set(levels) == {"trigger", "r1", "dl1", "dl2"}
    assert levels["trigger"]["event_type"] == EventType.SUBARRAY.value
    assert np.isclose(levels["trigger"]["time"], 1768964778.5)
    dl1 = levels["dl1"]
    assert dl1["camera_frame_hillas_intensity"] == 120.0 and np.isclose(dl1["camera_frame_hillas_psi"], 30.0)
    assert dl1["camera_frame_hillas_x"] == 0.05 and dl1["image"].shape == (1296,)
    dl2 = levels["dl2"]
    assert dl2["RF_tel_energy"] == 2.0 and dl2["Classifier_tel_prediction"] == 0.7
    assert dl2["Hillas_tel_impact_distance"] == 120.0
    stereo = messages[2][2]["dl2_stereo"]
    assert np.isclose(stereo["Hillas_alt"], np.rad2deg(1.2)) and stereo["Hillas_core_x"] == 10.0
    assert list(stereo["Hillas_telescopes"]) == [21, 22] and np.isclose(stereo["time"], 1768964778.5)
    # light events: no waveforms nor images; selection of the telescopes, no stereo
    light = event_messages(make_event(), full=False, tel_ids=[22], stereo=False)
    assert [topic for topic, _, _ in light] == ["tel22"]
    assert "r1" not in light[0][2] and "image" not in light[0][2]["dl1"]
    assert "camera_frame_hillas_intensity" in light[0][2]["dl1"]


def test_encode_decode():
    levels = {"r1": {"waveform": np.arange(6, dtype=np.int16).reshape(2, 3), "pixel_status": 3},
              "dl1": {"hillas_intensity": 5.5, "name": "x"}}
    frames = encode("tel21", levels, obs_id=1, event_id=2, tel_id=21, source="worker", count=4)
    assert frames[0] == b"tel21" and len(frames) == 3
    header, decoded = decode(frames)
    assert header["source"] == "worker" and header["count"] == 4
    assert np.array_equal(decoded["r1"]["waveform"], levels["r1"]["waveform"])
    assert decoded["r1"]["waveform"].dtype == np.int16
    assert decoded["dl1"] == {"hillas_intensity": 5.5, "name": "x", "obs_id": 1, "event_id": 2, "tel_id": 21}


def wait_for_publishers(receiver, duration=0.3):
    """
    Let the publishers connect: a SUB socket sends its subscriptions to a new publisher only when
    it is used (as the receiver of the viewer, which polls all the time)
    """
    end = time.monotonic() + duration
    while time.monotonic() < end:
        assert receiver.receive(timeout=0.02) is None


def receive_all(receiver, n, timeout=5.0):
    """The messages received (at most n, until timeout)"""
    messages, start = [], time.monotonic()
    while len(messages) < n and time.monotonic() - start < timeout:
        message = receiver.receive(timeout=0.1)
        if message is not None:
            messages.append(message)
    return messages


@pytest.mark.parametrize("n_publishers", [1, 2])
def test_publishers_to_receiver(n_publishers):
    """Several publishers (e.g. workers of the pipeline) connect to the receiver: each event once"""
    endpoint = "inproc://zmq-event-stream-test-%d" % n_publishers
    with ZMQEventReceiver(endpoint) as receiver:
        publishers = [ZMQEventPublisher(url=endpoint, full_event_period=0, source_name=f"worker{i}")
                      for i in range(n_publishers)]
        wait_for_publishers(receiver)
        for event_id in range(10):
            publishers[event_id % n_publishers](make_event(event_id), stereo=False)
        messages = receive_all(receiver, 20)
        for publisher in publishers:
            publisher.close()
    assert len(messages) == 20 and receiver.received == 20 and not receiver.lost
    received = sorted({(levels["dl1"]["event_id"], levels["dl1"]["tel_id"]) for _, levels in messages})
    assert received == [(event_id, tel_id) for event_id in range(10) for tel_id in (21, 22)]
    assert {header["source"] for header, _ in messages} == {f"worker{i}" for i in range(n_publishers)}


def test_topics_and_losses():
    endpoint = "inproc://zmq-event-stream-topics"
    with ZMQEventReceiver(endpoint, topics=["tel22"]) as receiver:
        publisher = ZMQEventPublisher(url=endpoint, full_event_period=3600)
        wait_for_publishers(receiver)
        publisher(make_event(1))
        publisher(make_event(2))
        messages = receive_all(receiver, 2)
        # only the topic tel22; the first event is full, the second light
        assert [levels["dl1"]["tel_id"] for _, levels in messages] == [22, 22]
        assert "r1" in messages[0][1] and "r1" not in messages[1][1]
        # the messages of the other topics are counted as lost (gaps of the counts of the publisher)
        assert receiver.lost == {publisher.source_name: 2}
        publisher.close()


def test_publisher_without_url():
    publisher = ZMQEventPublisher()
    assert not publisher.enabled
    publisher(make_event())  # nothing is sent


def test_process_publishes_the_events(tmp_path):
    """sst1mpipe-process --publish: the events written are received, with their DL1 parameters"""
    from ctapipe.io import read_table
    from conftest import TRANSITION_FILES, run_process

    tel_id = 21
    received, stop = [], threading.Event()
    with ZMQEventReceiver("tcp://127.0.0.1:*") as receiver:

        def receive():
            while not stop.is_set():
                message = receiver.receive(timeout=0.05)
                if message is not None:
                    received.append(message)

        thread = threading.Thread(target=receive, daemon=True)
        thread.start()
        output = tmp_path / "events.dl1.h5"
        run_process(TRANSITION_FILES[tel_id], output, f"--publish={receiver.endpoint}",
                    "--ZMQEventPublisher.full_event_period=0", "--ProcessorTool.allowed_sb_types=UNKNOWN")
        time.sleep(0.5)
        stop.set()
        thread.join()

    written = read_table(output, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
    by_event = {levels["dl1"]["event_id"]: levels for _, levels in received}
    # the first messages can be lost while the publisher connects: (almost) all the events are received
    assert len(by_event) >= 0.9 * len(written) and set(by_event) <= set(written["event_id"].tolist())
    for row in written:
        levels = by_event.get(int(row["event_id"]))
        if levels is None:
            continue
        assert levels["dl1"]["tel_id"] == tel_id
        np.testing.assert_allclose(levels["dl1"]["camera_frame_hillas_intensity"],
                                   row["camera_frame_hillas_intensity"])
        assert levels["dl1"]["image"].shape == (1296,) and "r0" in levels
        assert {"trigger", "r0", "r1", "dl1"} <= set(levels)
