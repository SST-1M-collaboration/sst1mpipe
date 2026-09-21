from typing import Dict, Generator
import zmq

from ctapipe.io import EventSource
from ctapipe.io.datalevels import DataLevel
from ctapipe.instrument import SubarrayDescription
from ctapipe.containers import SchedulingBlockContainer, ObservationBlockContainer, ArrayEventContainer
from protozfits import DL0v1_Telescope_pb2, CoreMessages_pb2, any_array_to_numpy, R1v1_pb2

from io.containers import SST1MArrayEventContainer

class ZMQEventSource(EventSource):


    def __init__(self, input_url, config=None, parent=None, **kwargs):

        super().__init__(input_url=input_url, config=config, parent=parent, **kwargs)
        context = zmq.Context()
        self.socket = context.socket(zmq.SUB)
        self.socket.connect(self.input_url)
        self.socket.subscribe(b"")

    @property
    def is_stream(self):
        return True

    @property
    def subarray(self) -> SubarrayDescription:
        """
        Obtain the subarray from the EventSource

        Returns
        -------
        ctapipe.instrument.SubarrayDecription

        """
        return None # TODO pass via config

    @property
    def observation_blocks(self) -> Dict[int, ObservationBlockContainer]:
        """
        Obtain the ObservationConfigurations from the EventSource, indexed by obs_id
        """
        raise NotImplementedError

    @property
    def scheduling_blocks(self) -> Dict[int, SchedulingBlockContainer]:
        """
        Obtain the ObservationConfigurations from the EventSource, indexed by obs_id
        """
        raise NotImplementedError

    @property
    def is_simulation(self) -> bool:
        """
        Whether the currently opened file is simulated

        Returns
        -------
        bool

        """
        return False

    @property
    def datalevels(self):

        return (DataLevel.R1, DataLevel.DL0)


    def _generator(self) -> Generator[ArrayEventContainer, None, None]:

        yield from self._generate_events()

    def _generate_events(self, event=None):

        if event is None:
            event = SST1MArrayEventContainer()
        data = self.socket.recv()

        msg = CoreMessages_pb2.CTAMessage()
        msg.ParseFromString(data)

        msg_type = msg.payload_type
        payload = msg.payload_data

        if msg_type == CoreMessages_pb2.DL0_TELESCOPE_EVENT:
            dl0_event =  DL0v1_Telescope_pb2.Event()
            dl0_event.ParseFromString(payload)

            tel_id = dl0_event.tel_id
            event.dl0.tel[tel_id].waveform = any_array_to_numpy(dl0_event.waveform)

        elif msg_type == CoreMessages_pb2.DL0_TELESCOPE_CAMERA_CONFIG:
            dl0_config = DL0v1_Telescope_pb2.CameraConfiguration()
            dl0_config.ParseFromString(payload)

            print("CAMERA CONFIG")
            print(dl0_config)

        elif msg_type == CoreMessages_pb2.DL0_TELESCOPE_DATA_STREAM:
            dl0_stream = DL0v1_Telescope_pb2.DataStream()
            dl0_stream.ParseFromString(payload)

            print("DATA STREAM")
            print(dl0_stream)

        elif msg_type == R1v1_pb2.TELESCOPE_DATA_STREAM:

            r1_stream = R1v1_pb2.TelescopeDataStream()
            r1_stream.ParseFromString(payload)

            pass
        elif msg_type == R1v1_pb2.CAMERA_CONFIG:

            r1_config = R1v1_pb2.CameraConfiguration()
            r1_config.ParseFromString(payload)
            pass

        elif msg_type == R1v1_pb2.R1_EVENT:

            r1_event = R1v1_pb2.Event()
            r1_event.ParseFromString(payload)
            pass

        if msg_type == R1v1_pb2.R1_EVENT or msg_type == DL0v1_Telescope_pb2.DL0_TELESCOPE_EVENT:

            return event

        return self._generate_events(event=event)



    def close(self):

        self.socket.close()
