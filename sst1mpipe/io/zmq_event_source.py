from typing import Dict

import numpy as np
import zmq
from astropy.time import Time
from importlib.resources import files


from ctapipe.io import EventSource
from ctapipe.io.datalevels import DataLevel
from ctapipe.instrument import SubarrayDescription
from ctapipe.containers import SchedulingBlockContainer, ObservationBlockContainer, DL0Container, R1Container
from protozfits import DL0v1_Telescope_pb2, CoreMessages_pb2, any_array_to_numpy, R1v1_pb2
from ctapipe.core.traits import Unicode, Path

from sst1mpipe.io.containers import SST1MArrayEventContainer

def ctao_high_res_to_time(seconds, quarter_nanoseconds): # TODO import from ctapipe==0.24
    """Convert CTAO high resolution timestamp to astropy Time."""
    # unix_tai accepts two floats for maximum precision
    # we can just pass integral and fractional part
    fractional_seconds = quarter_nanoseconds * 0.25e-9
    return Time(
        seconds,
        fractional_seconds,
        format="unix_tai",
        # this is only for displaying iso timestamp, not any actual precision
        precision=9,
    )

def fill_DL0v1_Telescope_Event_to_DL0Container(payload: bytes, dl0: DL0Container) -> int:

    dl0_message = DL0v1_Telescope_pb2.Event()
    dl0_message.ParseFromString(payload)

    n_chan, n_pix, n_samples = (dl0_message.num_channels, dl0_message.num_pixels_survived, dl0_message.num_samples)

    tel_id = dl0_message.tel_id
    dl0.tel[tel_id].event_type = dl0_message.event_type
    dl0.tel[tel_id].event_time = ctao_high_res_to_time(dl0_message.event_time_s, dl0_message.event_time_qns) # TODO use ctapipe > 0.24 with ctapipe.time.ctao_high_res_to_time
    dl0.tel[tel_id].waveform = any_array_to_numpy(dl0_message.waveform).reshape((n_chan, n_pix, n_samples)) - any_array_to_numpy(dl0_message.pedestal_intensity).reshape((n_chan, n_pix,))[..., np.newaxis]
    dl0.tel[tel_id].pixel_status = any_array_to_numpy(dl0_message.pixel_status)
    dl0.tel[tel_id].first_cell_id = any_array_to_numpy(dl0_message.first_cell_id)
    dl0.tel[tel_id].calibration_monitoring_id = dl0_message.calibration_monitoring_id

    return dl0_message.event_id

def fill_R1v1_Event_to_R1Container(payload: bytes, r1: R1Container) -> int:

    r1_message = R1v1_pb2.Event()
    r1_message.ParseFromString(payload)

    n_chan, n_pix, n_samples = (r1_message.num_channels, r1_message.num_pixels, r1_message.num_samples)

    tel_id = r1_message.tel_id
    r1.tel[tel_id].event_type = r1_message.event_type
    r1.tel[tel_id].event_time = ctao_high_res_to_time(r1_message.event_time_s, r1_message.event_time_qns) # TODO use ctapipe > 0.24 with ctapipe.time.ctao_high_res_to_time
    r1.tel[tel_id].waveform = any_array_to_numpy(r1_message.waveform).reshape((n_chan, n_pix, n_samples))
    r1.tel[tel_id].pedestal_intensity = any_array_to_numpy(r1_message.pedestal_intensity).reshape((n_chan, n_pix))
    r1.tel[tel_id].pixel_status = any_array_to_numpy(r1_message.pixel_status)
    r1.tel[tel_id].first_cell_id = any_array_to_numpy(r1_message.first_cell_id)
    r1.tel[tel_id].module_hires_local_clock_counter = any_array_to_numpy(r1_message.module_hires_local_clock_counter)
    r1.tel[tel_id].calibration_monitoring_id = r1_message.calibration_monitoring_id

    return r1_message.event_id, r1_message.local_run_id

class ZMQEventSource(EventSource):

    input_url = Unicode(info_text="URL of the input stream",
                        help="TCP and port address for the input ZMQ stream. Example `tcp://192.168.1.1:1986` ")

    subarray_file = Path(help="Path to the file containing the subarray-description.",
                         default_value=files('sst1mpipe.data').joinpath('sst1m_array.h5')).tag(config=True)

    def __init__(self, input_url, config=None, parent=None, **kwargs):

        super().__init__(input_url=input_url, config=config, parent=parent, **kwargs)
        context = zmq.Context()
        self.socket = context.socket(zmq.PULL)
        self.socket.connect(self.input_url)
        self._subarray = SubarrayDescription.from_hdf(self.subarray_file)
        self._scheduling_blocks = {tel_id: SchedulingBlockContainer() for tel_id in self.subarray.tel_ids}
        self._observations_blocks = {tel_id: ObservationBlockContainer() for tel_id in self.subarray.tel_ids}

    @staticmethod
    def is_compatible(file_path: str) -> bool:

        context = zmq.Context()
        socket = context.socket(zmq.PULL)

        try:
            socket.connect(file_path)
            return True
        except zmq.ZMQError:
            return False
        finally:
            try:
                socket.close(linger=0)
            except Exception:
                pass
            context.term()

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
        return self._subarray

    @property
    def observation_blocks(self) -> Dict[int, ObservationBlockContainer]:
        """
        Obtain the ObservationConfigurations from the EventSource, indexed by obs_id
        """
        UserWarning("Observations blocks is not yet implemented returns default empty containers")
        return self._observations_blocks

    @property
    def scheduling_blocks(self) -> Dict[int, SchedulingBlockContainer]:
        """
        Obtain the ObservationConfigurations from the EventSource, indexed by obs_id
        """
        UserWarning("Scheduling blocks is not yet implemented returns default empty containers")
        return self._scheduling_blocks

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


    def _generator(self):

        yield from self._generate_events()

    def _generate_events(self):

        count = 0

        while True:

            event = SST1MArrayEventContainer()
            data = self.socket.recv()

            msg = CoreMessages_pb2.CTAMessage()
            msg.ParseFromString(data)

            msg_types = msg.payload_type
            payloads = msg.payload_data

            for msg_type, payload in zip(msg_types, payloads, strict=True):
                if msg_type == CoreMessages_pb2.DL0_TELESCOPE_EVENT:

                    event_id = fill_DL0v1_Telescope_Event_to_DL0Container(payload, event.dl0)
                    event.index.event_id = event_id

                elif msg_type == CoreMessages_pb2.DL0_TELESCOPE_CAMERA_CONFIG:
                    dl0_config = DL0v1_Telescope_pb2.CameraConfiguration()
                    dl0_config.ParseFromString(payload)

                    raise NotImplementedError

                elif msg_type == CoreMessages_pb2.DL0_TELESCOPE_DATA_STREAM:
                    dl0_stream = DL0v1_Telescope_pb2.DataStream()
                    dl0_stream.ParseFromString(payload)

                    raise NotImplementedError

                elif msg_type == CoreMessages_pb2.TELESCOPE_DATA_STREAM:

                    r1_stream = R1v1_pb2.TelescopeDataStream()
                    r1_stream.ParseFromString(payload)

                    raise NotImplementedError

                elif msg_type == CoreMessages_pb2.CAMERA_CONFIG:

                    r1_config = R1v1_pb2.CameraConfiguration()
                    r1_config.ParseFromString(payload)

                    raise NotImplementedError

                elif msg_type == CoreMessages_pb2.R1_EVENT:

                    event_id, local_id = fill_R1v1_Event_to_R1Container(payload, event.r1) # TODO use local_id
                    event.index.event_id = event_id


                if msg_type == CoreMessages_pb2.R1_EVENT or msg_type == CoreMessages_pb2.DL0_TELESCOPE_EVENT:

                    event.count = count
                    count += 1
                    yield event



    def close(self):

        self.socket.close()
