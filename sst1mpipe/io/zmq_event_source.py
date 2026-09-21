from typing import Dict, Generator
import zmq
import ipaddress

from ctapipe.io import EventSource
from ctapipe.io.datalevels import DataLevel
from ctapipe.instrument import SubarrayDescription
from ctapipe.containers import SchedulingBlockContainer, ObservationBlockContainer, ArrayEventContainer, \
    DL0Container
from protozfits import DL0v1_Telescope_pb2, CoreMessages_pb2, any_array_to_numpy, R1v1_pb2
from ctapipe.core.traits import Unicode

from sst1mpipe.io.containers import SST1MArrayEventContainer

def fill_DL0v1_Telescope_Event_to_DL0Container(payload: bytes, dl0: DL0Container) -> DL0Container:

    dl0_event = DL0v1_Telescope_pb2.Event()
    dl0_event.ParseFromString(payload)

    tel_id = dl0_event.tel_id
    dl0.tel[tel_id].waveform = any_array_to_numpy(dl0_event.waveform)
    


class ZMQEventSource(EventSource):

    input_url = Unicode(info_text="URL of the input stream",
                        help="TCP and port address for the input ZMQ stream. Example `tcp://192.168.1.1:1986` ")

    def __init__(self, input_url, config=None, parent=None, **kwargs):

        super().__init__(input_url=input_url, config=config, parent=parent, **kwargs)
        context = zmq.Context()
        self.socket = context.socket(zmq.PULL)
        self.socket.connect(self.input_url)

    def is_compatible(self, file_path: str) -> bool:

        try:
            protocol, host, port = file_path.split(":", 2)
        except ValueError:
            return False

        if protocol != "tcp":

            return False

        try:
            ipaddress.ip_address(host)
        except ValueError:
            return False

        try:
            port = int(port)
        except ValueError:
            return False

        if not 1 <= port <= 65535:
            return False

        return True

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

        msg_types = msg.payload_type
        payloads = msg.payload_data

        for msg_type, payload in zip(msg_types, payloads):


            if msg_type == CoreMessages_pb2.DL0_TELESCOPE_EVENT:
                
                fill_DL0v1_Telescope_Event_to_DL0Container(payload, event.dl0)

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

            elif msg_type == CoreMessages_pb2.TELESCOPE_DATA_STREAM:

                r1_stream = R1v1_pb2.TelescopeDataStream()
                r1_stream.ParseFromString(payload)

                pass
            elif msg_type == CoreMessages_pb2.CAMERA_CONFIG:

                r1_config = R1v1_pb2.CameraConfiguration()
                r1_config.ParseFromString(payload)
                pass

            elif msg_type == CoreMessages_pb2.R1_EVENT:

                r1_event = R1v1_pb2.Event()
                r1_event.ParseFromString(payload)
                pass

            if msg_type == CoreMessages_pb2.R1_EVENT or msg_type == CoreMessages_pb2.DL0_TELESCOPE_EVENT:

                return event

            return self._generate_events(event=event)



    def close(self):

        self.socket.close()
