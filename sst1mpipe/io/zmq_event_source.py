"""
Event source of the telescope events sent in a ZMQ stream (protocol buffers, e.g. by the
camera server for the real time analysis): DigiCam R0 events (``DataModel.CameraEvent``,
as in the zfits files) and CTAO R1 and DL0 events.
"""
from dataclasses import dataclass
from functools import partial

import numpy as np
import zmq
from ctapipe.containers import (
    EventType,
    ObservationBlockContainer,
    SchedulingBlockContainer,
    SchedulingBlockType,
)
from ctapipe.core.traits import Dict, List, Path, Undefined, Unicode, UseEnum
from ctapipe.instrument import SubarrayDescription
from ctapipe.io import EventSource
from ctapipe.io.datalevels import DataLevel
from ctapipe.time import ctao_high_res_to_time
from protozfits import (
    CoreMessages_pb2,
    DL0v1_Telescope_pb2,
    ProtoDataModel_pb2,
    R1v1_pb2,
    any_array_to_numpy,
    make_namedtuple,
)

from sst1mpipe.io.containers import SST1MArrayEventContainer, StreamCameraConfigContainer
from sst1mpipe.io.sst1m_event_source import fill_r0_event, parse_swapped_modules, pixel_order
from sst1mpipe.resources import SUBARRAY_FILE

__all__ = ["ZMQEventSource", "DataStream"]

CAMERA_CONFIG_FIELDS = (
    "tel_id", "local_run_id", "config_time_s", "camera_config_id", "num_modules", "num_pixels",
    "num_channels", "num_samples_nominal", "num_samples_long", "num_samples_removed_start",
    "num_samples_removed_end", "data_model_version", "calibration_service_id", "calibration_algorithm_id",
)


@dataclass
class DataStream:
    """Data stream of a telescope: its scheduling and observation blocks and the waveform encoding"""

    tel_id: int
    sb_id: int
    obs_id: int
    waveform_scale: float = 1.0
    waveform_offset: float = 0.0


def event_type(value):
    """ctapipe EventType of the event type of a message, UNKNOWN if it is not a known type"""
    try:
        return EventType(value)
    except ValueError:
        return EventType.UNKNOWN


def decode_waveform(stored, data_stream):
    """
    Waveform from its stored values: ``stored / waveform_scale - waveform_offset``
    (CTAO R1/DL0 data model). Without data stream, the stored values are used.
    """
    waveform = stored.astype(np.float32)
    if data_stream is None:
        return waveform
    return waveform / data_stream.waveform_scale - data_stream.waveform_offset


def camera_config_container(message, data_level):
    """StreamCameraConfigContainer of a R1 or DL0 CameraConfiguration message"""
    config = StreamCameraConfigContainer(
        data_level=data_level,
        pixel_id_map=any_array_to_numpy(message.pixel_id_map),
        module_id_map=any_array_to_numpy(message.module_id_map),
        **{name: getattr(message, name) for name in CAMERA_CONFIG_FIELDS},
    )
    if data_level == "DL0":
        config.sampling_frequency = message.sampling_frequency
    return config


class ZMQEventSource(EventSource):
    """
    Read the telescope events of a ZMQ stream (PULL socket).

    The messages of the stream (``CTAMessage``) are:

    - the DigiCam camera events (``CAMERA_EVENT``, ``DataModel.CameraEvent`` as in the zfits
      files), filling ``event.r0`` and the trigger as the `SST1MEventSource`. The event id is
      the SWAT array event id if the first event has one, the camera event number otherwise.
    - the data stream of a telescope (``TELESCOPE_DATA_STREAM``, ``DL0_TELESCOPE_DATA_STREAM``):
      its scheduling and observation blocks and the encoding of the waveforms
      (``stored / waveform_scale - waveform_offset``)
    - the camera configuration (``CAMERA_CONFIG``, ``DL0_TELESCOPE_CAMERA_CONFIG``),
      available in ``camera_config``
    - the events (``R1_EVENT``, ``DL0_TELESCOPE_EVENT``), filling ``event.r1`` or ``event.dl0``,
      the event index and the trigger
    - ``END_OF_STREAM``, which ends the reading.

    The scheduling and observation blocks are only known once the data stream message is
    received; the type of the scheduling block is not sent and is ``scheduling_block_type``.
    """

    input_url = Unicode(
        info_text="URL of the input stream",
        help="TCP and port address for the input ZMQ stream. Example `tcp://192.168.1.1:1986` ",
    ).tag(config=True)

    subarray_file = Path(
        help="Path to the file containing the subarray-description.",
        default_value=SUBARRAY_FILE,
    ).tag(config=True)

    scheduling_block_type = UseEnum(
        SchedulingBlockType,
        default_value=SchedulingBlockType.OBSERVATION,
        help="Type of the scheduling blocks of the stream, which is not sent in the stream",
    ).tag(config=True)

    swapped_modules = List(
        trait=Dict(),
        default_value=[],
        help=(
            "Wrongly connected modules of the DigiCam camera events, see"
            " SST1MEventSource.swapped_modules"
        ),
    ).tag(config=True)

    def __init__(self, input_url=Undefined, config=None, parent=None, **kwargs):

        super().__init__(input_url=input_url, config=config, parent=parent, **kwargs)
        self._context = zmq.Context()
        self.socket = self._context.socket(zmq.PULL)
        self.socket.connect(self.input_url)
        self.log.info("Connect to socket %s", self.input_url)
        self._subarray = SubarrayDescription.from_hdf(self.subarray_file)
        self._scheduling_blocks = {}
        self._observation_blocks = {}
        self._data_streams = {}
        self._camera_configs = {}
        self._end_of_stream = False
        self._pixel_order = partial(pixel_order, parse_swapped_modules(self.swapped_modules))
        self._swat_event_ids = None

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
        return self._subarray

    @property
    def observation_blocks(self) -> dict[int, ObservationBlockContainer]:
        """Observation blocks of the data streams received, by obs_id"""
        return self._observation_blocks

    @property
    def scheduling_blocks(self) -> dict[int, SchedulingBlockContainer]:
        """Scheduling blocks of the data streams received, by sb_id"""
        return self._scheduling_blocks

    @property
    def data_streams(self) -> dict[int, DataStream]:
        """Data streams received, by telescope id"""
        return self._data_streams

    @property
    def camera_config(self) -> dict[int, StreamCameraConfigContainer]:
        """Camera configurations received, by telescope id"""
        return self._camera_configs

    @property
    def is_simulation(self) -> bool:
        return False

    @property
    def datalevels(self):
        return (DataLevel.R0, DataLevel.R1, DataLevel.DL0)

    def _generator(self):
        count = 0
        while not self._end_of_stream:
            message = CoreMessages_pb2.CTAMessage()
            message.ParseFromString(self.socket.recv())
            for msg_type, payload in zip(message.payload_type, message.payload_data, strict=True):
                event = self._read_payload(msg_type, payload)
                if event is not None:
                    event.count = count
                    count += 1
                    yield event
                if self._end_of_stream:
                    break

    def _read_payload(self, msg_type, payload):
        """Read a payload of a message: the event for an event payload, None otherwise"""
        if msg_type == CoreMessages_pb2.CAMERA_EVENT:
            return self._camera_event(payload)
        if msg_type == CoreMessages_pb2.R1_EVENT:
            return self._r1_event(payload)
        if msg_type == CoreMessages_pb2.DL0_TELESCOPE_EVENT:
            return self._dl0_event(payload)
        if msg_type == CoreMessages_pb2.TELESCOPE_DATA_STREAM:
            stream = R1v1_pb2.TelescopeDataStream()
            stream.ParseFromString(payload)
            self._add_data_stream(stream, "R1")
        elif msg_type == CoreMessages_pb2.DL0_TELESCOPE_DATA_STREAM:
            stream = DL0v1_Telescope_pb2.DataStream()
            stream.ParseFromString(payload)
            self._add_data_stream(stream, "DL0")
        elif msg_type == CoreMessages_pb2.CAMERA_CONFIG:
            config = R1v1_pb2.CameraConfiguration()
            config.ParseFromString(payload)
            self._add_camera_config(config, "R1")
        elif msg_type == CoreMessages_pb2.DL0_TELESCOPE_CAMERA_CONFIG:
            config = DL0v1_Telescope_pb2.CameraConfiguration()
            config.ParseFromString(payload)
            self._add_camera_config(config, "DL0")
        elif msg_type == CoreMessages_pb2.END_OF_STREAM:
            self.log.info("End of stream received")
            self._end_of_stream = True
        else:
            self.log.debug("Message of type %d ignored", msg_type)
        return None

    def _add_data_stream(self, stream, data_level):
        tel_id = stream.tel_id
        data_stream = DataStream(
            tel_id=tel_id, sb_id=stream.sb_id, obs_id=stream.obs_id,
            waveform_scale=stream.waveform_scale or 1.0, waveform_offset=stream.waveform_offset,
        )
        self._data_streams[tel_id] = data_stream
        producer_id = f"SST1M-{tel_id}"
        self._scheduling_blocks[stream.sb_id] = SchedulingBlockContainer(
            sb_id=np.uint64(stream.sb_id), sb_type=self.scheduling_block_type, producer_id=producer_id,
        )
        self._observation_blocks[stream.obs_id] = ObservationBlockContainer(
            obs_id=np.uint64(stream.obs_id), sb_id=np.uint64(stream.sb_id), producer_id=producer_id,
        )
        self.log.info(
            "%s data stream of telescope %d: sb_id %d, obs_id %d, waveform scale %.2f, offset %.2f",
            data_level, tel_id, stream.sb_id, stream.obs_id, data_stream.waveform_scale, data_stream.waveform_offset,
        )

    def _add_camera_config(self, message, data_level):
        config = camera_config_container(message, data_level)
        self._camera_configs[config.tel_id] = config
        self.log.info(
            "%s camera configuration of telescope %d, local run id %d: data shape (%d, %d, %d)",
            data_level, config.tel_id, config.local_run_id,
            config.num_channels, config.num_pixels, config.num_samples_nominal,
        )

    def _fill_index_and_trigger(self, event, tel_id, event_id, obs_id, camera):
        event.index.event_id = event_id
        event.index.obs_id = obs_id
        event.trigger.event_type = camera.event_type
        event.trigger.time = camera.event_time
        event.trigger.tel[tel_id].time = camera.event_time
        event.trigger.tels_with_trigger = [tel_id]

    def _camera_event(self, payload):
        """R0 event of a DigiCam camera event (DataModel.CameraEvent)"""
        message = ProtoDataModel_pb2.CameraEvent()
        message.ParseFromString(payload)
        camera_event = make_namedtuple(message)

        event = SST1MArrayEventContainer()
        event.r0.meta = dict(is_simulation=False)
        tel_id = fill_r0_event(event, camera_event, self._pixel_order)

        # the event ids of the stream are the SWAT array event ids if the first event has one
        if self._swat_event_ids is None:
            self._swat_event_ids = camera_event.arrayEvtNum != 0
            self.log.info(
                "Event id of the camera events: %s",
                "SWAT array event id" if self._swat_event_ids else "camera event number",
            )
        event.index.event_id = camera_event.arrayEvtNum if self._swat_event_ids else camera_event.eventNumber
        data_stream = self._data_streams.get(tel_id)
        event.index.obs_id = data_stream.obs_id if data_stream is not None else 0
        return event

    def _r1_event(self, payload):
        message = R1v1_pb2.Event()
        message.ParseFromString(payload)
        tel_id = message.tel_id
        n_channels, n_pixels, n_samples = message.num_channels, message.num_pixels, message.num_samples
        data_stream = self._data_streams.get(tel_id)

        event = SST1MArrayEventContainer()
        r1 = event.r1.tel[tel_id]
        r1.event_type = event_type(message.event_type)
        r1.event_time = ctao_high_res_to_time(message.event_time_s, message.event_time_qns)
        stored = any_array_to_numpy(message.waveform).reshape((n_channels, n_pixels, n_samples))
        r1.waveform = decode_waveform(stored, data_stream)
        r1.pedestal_intensity = any_array_to_numpy(message.pedestal_intensity).reshape((n_channels, n_pixels))
        r1.pixel_status = any_array_to_numpy(message.pixel_status)
        r1.first_cell_id = any_array_to_numpy(message.first_cell_id)
        r1.module_hires_local_clock_counter = any_array_to_numpy(message.module_hires_local_clock_counter)
        r1.calibration_monitoring_id = message.calibration_monitoring_id
        r1.selected_gain_channel = np.zeros(n_pixels, dtype=int)

        # obs_id of the data stream, or the local run id of the camera before the data stream
        obs_id = data_stream.obs_id if data_stream is not None else message.local_run_id
        self._fill_index_and_trigger(event, tel_id, message.event_id, obs_id, r1)
        return event

    def _dl0_event(self, payload):
        message = DL0v1_Telescope_pb2.Event()
        message.ParseFromString(payload)
        tel_id = message.tel_id
        n_channels, n_pixels, n_samples = message.num_channels, message.num_pixels_survived, message.num_samples
        data_stream = self._data_streams.get(tel_id)

        event = SST1MArrayEventContainer()
        dl0 = event.dl0.tel[tel_id]
        dl0.event_type = event_type(message.event_type)
        dl0.event_time = ctao_high_res_to_time(message.event_time_s, message.event_time_qns)
        stored = any_array_to_numpy(message.waveform).reshape((n_channels, n_pixels, n_samples))
        pedestal = any_array_to_numpy(message.pedestal_intensity).reshape((n_channels, n_pixels))
        dl0.waveform = decode_waveform(stored, data_stream) - pedestal[..., np.newaxis]
        dl0.pixel_status = any_array_to_numpy(message.pixel_status)
        dl0.first_cell_id = any_array_to_numpy(message.first_cell_id)
        dl0.calibration_monitoring_id = message.calibration_monitoring_id
        dl0.selected_gain_channel = np.zeros(n_pixels, dtype=int)

        obs_id = data_stream.obs_id if data_stream is not None else 0
        self._fill_index_and_trigger(event, tel_id, message.event_id, obs_id, dl0)
        return event

    def close(self):
        self.socket.close(linger=0)
        self._context.term()
        self.log.info("Closed ZMQ connection to %s", self.input_url)
