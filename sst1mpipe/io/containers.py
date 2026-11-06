"""
Containers of the SST-1M data. The events are ctapipe `~ctapipe.containers.ArrayEventContainer`
(`SST1MArrayEventContainer`) with the DigiCam specific raw data (R0) and the statistics of
the ADC samples of the pedestal events in the camera monitoring. The other data levels follow ctapipe.
The PyTables descriptions are those of the tables written by sst1mpipe in the DL1/DL2 files.
"""
from enum import Flag, IntEnum
from functools import partial

import numpy as np
from astropy.time import TimeDelta
from ctapipe.containers import (
    ArrayEventContainer,
    CameraMonitoringContainer,
    ChunkStatisticsContainer,
    EventType,
    MonitoringContainer,
    ObservationBlockContainer,
    PixelStatisticsContainer,
    R0CameraContainer,
    R0Container,
    TelescopeMonitoringContainer,
)
from ctapipe.core import Container, Field, Map

from tables import (
    BoolCol,
    Float64Col,
    Int64Col,
    IsDescription,
    StringCol,
)

__all__ = [
    "CameraEventType",
    "DigicamCardType",
    "PEDESTAL_EVENT_TYPES",
    "SST1MR0CameraContainer",
    "SST1MR0Container",
    "SST1MPixelStatisticsContainer",
    "SST1MCameraMonitoringContainer",
    "SST1MTelescopeMonitoringContainer",
    "SST1MMonitoringContainer",
    "SST1MArrayEventContainer",
    "SST1MObservationBlockContainer",
    "DigicamConfigContainer",
    "StreamCameraConfigContainer",
    "DL1_info",
    "DL2_info",
]


class CameraEventType(Flag):
    """
    Event type of a DigiCam camera event: a bit mask of the trigger conditions of the event.

    Each trigger condition is one bit and the event type is the sum of the bits of the event:

    =====  ===========  ================  ==================================================
    Bit    Value        Member            Meaning
    =====  ===========  ================  ==================================================
    0      ``0x1``      ``PATCH7``        algorithm 0 trigger (patch7)
    1      ``0x2``      ``PATCH19``       algorithm 1 trigger (patch19)
    2      ``0x4``      ``MUON_TRIGGER``  algorithm 2 trigger (muon)
    3      ``0x8``      ``INTERNAL``      internal trigger (the pedestal events)
    4      ``0x10``     ``EXTMSTR``       external trigger, on the master only
    5, 6   ``0x20``,    ``BIT5``,         unused
           ``0x40``     ``BIT6``
    7      ``0x80``     ``CONTINUOUS``    continuous readout marker
    16     ``0x10000``  ``MUON_DETECT``   muon detected by the camera server
    17     ``0x20000``  ``HILLAS``        Hillas parameters computed by the camera server
    =====  ===========  ================  ==================================================

    E.g. ``5 = 0x1 + 0x4`` is an event triggered by both the patch7 and the muon triggers,
    ``PATCH7|MUON_TRIGGER``.

    The readers (``SST1MEventSource``, ``ZMQEventSource``) store it in
    ``event.r0.tel[tel_id].event_type``: the container converts the integer of the file or
    stream to a ``CameraEventType``.

    Notes
    -----
    - ``==`` compares the whole mask: ``CameraEventType(5) == CameraEventType.PATCH7`` is False.
      To test if an event has a trigger condition, whatever the other bits, use ``in`` (or ``&``).
    - ``UNKNOWN`` (0) is in every event type, as the empty mask of a ``Flag``: test
      ``not event_type`` for an event without any bit.
    - A bit which is not a member (bits 8 to 15, above 17) raises a ``ValueError``.
    - The ctapipe event type of the event (``event.trigger.event_type``) is a pedestal event for
      an event type equal to ``INTERNAL``, ``SUBARRAY`` otherwise, see
      ``sst1mpipe.io.sst1m_event_source.camera_event_type_to_event_type``.

    Examples
    --------
    Test the trigger conditions of an event, e.g. of type 5:

    >>> from sst1mpipe.io.containers import CameraEventType
    >>> event_type = CameraEventType(5)  # e.g. event.r0.tel[tel_id].event_type
    >>> event_type
    <CameraEventType.PATCH7|MUON_TRIGGER: 5>
    >>> CameraEventType.PATCH7 in event_type, CameraEventType.MUON_TRIGGER in event_type
    (True, True)
    >>> CameraEventType.PATCH19 in event_type
    False
    >>> event_type == CameraEventType.PATCH7
    False

    Any of several trigger conditions:

    >>> bool(event_type & (CameraEventType.PATCH7 | CameraEventType.PATCH19))
    True

    Bits, name and value of the event type:

    >>> list(event_type)
    [<CameraEventType.PATCH7: 1>, <CameraEventType.MUON_TRIGGER: 4>]
    >>> event_type.name, event_type.value
    ('PATCH7|MUON_TRIGGER', 5)
    >>> CameraEventType.PATCH7 | CameraEventType.MUON_TRIGGER == event_type
    True

    The integers of a table (e.g. the ``event_type`` column of ``sst1mpipe-trigger``) are tested
    with the values of the bits:

    >>> import numpy as np
    >>> event_types = np.array([1, 5, 8, 4])
    >>> (event_types & CameraEventType.MUON_TRIGGER.value) != 0
    array([False,  True, False,  True])
    """

    UNKNOWN = 0x0
    PATCH7 = 0x1  # algorithm 0 trigger - PATCH7
    PATCH19 = 0x2  # algorithm 1 trigger - PATCH19
    MUON_TRIGGER = 0x4  # algorithm 2 trigger - MUON
    INTERNAL = 0x8  # internal or external trigger
    EXTMSTR = 0x10  # unused (0) / external (on master only)
    BIT5 = 0x20  # unused (0)
    BIT6 = 0x40  # unused (0)
    CONTINUOUS = 0x80  # continuous readout marker
    MUON_DETECT = 0x10000  # camera server detected muon
    HILLAS = 0x20000  # camera server computed Hillas parametrs


class DigicamCardType(IntEnum):
    """Type of a DigiCam board (``DigicamCardType`` enum of the ``DigicamConfig`` protobuf message)"""

    TRIGGER = 0
    FADC_ANALOG = 1
    FADC_INTERSIL = 2


# ctapipe event types of the SST-1M pedestal events: SKY_PEDESTAL, DARK_PEDESTAL in the dark runs
# and ELECTRONIC_PEDESTAL in the bias runs
PEDESTAL_EVENT_TYPES = (EventType.SKY_PEDESTAL, EventType.DARK_PEDESTAL, EventType.ELECTRONIC_PEDESTAL)


class SST1MR0CameraContainer(R0CameraContainer):
    """
    Raw data of a single SST-1M telescope: the ctapipe
    `~ctapipe.containers.R0CameraContainer` (``waveform`` of shape
    (n_channels, n_pixels, n_samples), n_channels = 1 for DigiCam)
    with the DigiCam specific information.

    The static information of the camera (geometry, number of pixels and samples)
    is in the `~ctapipe.instrument.SubarrayDescription` of the event source,
    the trigger cluster and patch matrices in `sst1mpipe.instrument.camera.DigiCam`.
    """

    pixel_flags = Field(None, "numpy array containing pixel flags (n_pixels)")
    pedestal = Field(None, "baseline computed by DigiCam from 1024 pre-samples, in ADC (n_pixels), 0 for the pixels without data")
    camera_event_number = Field(None, "event number within the first trigger of operation")
    event_time = Field(None, "White Rabbit timestamp of the event (ns, TAI)")
    gps_time = Field(None, "timestamp from a precise external clock (ns)")
    _event_type = Field(None, "camera event type")
    trigger_input_traces = Field(None, "trigger patch traces (n_patches, n_samples), patch_sw_id order")
    trigger_output_patch7 = Field(None, "trigger 7 patch cluster traces (n_clusters, n_samples), central patch_sw_id order")
    trigger_output_patch19 = Field(None, "trigger 19 patch cluster traces (n_clusters, n_samples), central patch_sw_id order")
    trigger_output_muon = Field(None, "trigger muon cluster traces (n_clusters, n_samples), central patch_sw_id order")

    @property
    def event_type(self):
        return self._event_type

    @event_type.setter
    def event_type(self, value):
        self._event_type = CameraEventType(value)


class SST1MR0Container(R0Container):
    """
    Raw data of the SST-1M telescopes
    """

    tel = Field(
        default_factory=partial(Map, SST1MR0CameraContainer),
        description="map of tel_id to SST1MR0CameraContainer",
    )


class SST1MPixelStatisticsContainer(PixelStatisticsContainer):
    """
    ctapipe pixel statistics with the statistics of the ADC samples of the pedestal events.
    The statistics of the pedestal events are computed in chunks of events
    (`~ctapipe.containers.ChunkStatisticsContainer`) by sst1mpipe-monitoring, and filled by
    `sst1mpipe.io.monitoring_source.SST1MMonitoringSource`.
    """

    pedestal_image = Field(
        default_factory=ChunkStatisticsContainer,
        description=(
            "Statistics of the calibrated images (p.e.) of the pedestal events, shape (n_pixels)."
            " Used by the NSBImageCleaner"
        ),
    )
    pedestal_waveform = Field(
        default_factory=ChunkStatisticsContainer,
        description=(
            "Statistics of the ADC samples of the pedestal events, shape (n_channels, n_pixels):"
            " mean and median over the events of the mean of the samples (baseline), and mean over"
            " the events of the std of the samples (used for the voltage drop correction and the"
            " dead pixels)"
        ),
    )


class SST1MCameraMonitoringContainer(CameraMonitoringContainer):
    """ctapipe camera monitoring with the SST-1M pixel statistics"""

    pixel_statistics = Field(
        default_factory=SST1MPixelStatisticsContainer,
        description="Pixel statistics of the pedestal events",
    )


class SST1MTelescopeMonitoringContainer(TelescopeMonitoringContainer):
    """ctapipe telescope monitoring with the SST-1M camera monitoring"""

    camera = Field(
        default_factory=SST1MCameraMonitoringContainer,
        description="Monitoring data of the camera",
    )


class SST1MMonitoringContainer(MonitoringContainer):
    tel = Field(
        default_factory=partial(Map, SST1MTelescopeMonitoringContainer),
        description="map of tel_id to SST1MTelescopeMonitoringContainer",
    )


class DigicamConfigContainer(Container):
    """Configuration of the DigiCam"""

    default_prefix = "digicam"

    protocol_vers = Field(None, "Protocol version of each board (n_boards)")
    sn = Field(None, "Serial number of each board (n_boards), 0 for the empty board slots")
    hv = Field(None, "Hardware version of each board (n_boards)")
    gateware_rev = Field(None, "Gateware revision of each board (n_boards)")
    gateware_vers = Field(None, "Gateware version of each board (n_boards)")
    gateware_code = Field(None, "Gateware code of each board (n_boards)")
    gateware_card_type = Field(None, "Gateware card type of each board, `DigicamCardType` (n_boards)")
    firmware_rev = Field(None, "Firmware revision of each board (n_boards)")
    firmware_vers = Field(None, "Firmware version of each board (n_boards)")
    firmware_code = Field(None, "Firmware code of each board (n_boards)")
    firmware_card_type = Field(None, "Firmware card type of each board, `DigicamCardType` (n_boards)")
    operation_id = Field(-1, "Operation id (reserved)")
    operation_data = Field(-1, "Operation data (reserved)")
    digicam_time_sec = Field(-1, "DigiCam internal clock, seconds")
    digicam_time_nanosec = Field(-1, "DigiCam internal clock, nanoseconds")

    @property
    def digicam_time(self):
        """
        DigiCam internal clock as an astropy `~astropy.time.TimeDelta`, with ns precision:
        a duration since the origin of the clock, without absolute reference
        """
        return TimeDelta(self.digicam_time_sec, self.digicam_time_nanosec / 1e9, format="sec")


class StreamCameraConfigContainer(Container):
    """
    Configuration of the camera of a telescope, sent in a ZMQ stream
    (R1 ``CameraConfiguration`` or DL0 ``Telescope.CameraConfiguration`` message)
    """

    default_prefix = "camera_config"

    data_level = Field("", "Data level of the stream, R1 or DL0")
    tel_id = Field(-1, "Telescope id")
    local_run_id = Field(-1, "Local run id of the camera")
    config_time_s = Field(np.nan, "Time of the configuration (s)")
    camera_config_id = Field(-1, "Id of the camera configuration")
    pixel_id_map = Field(None, "Pixel id of each pixel of the data (n_pixels)")
    module_id_map = Field(None, "Module id of each module of the data (n_modules)")
    num_modules = Field(-1, "Number of modules")
    num_pixels = Field(-1, "Number of pixels")
    num_channels = Field(-1, "Number of gain channels")
    num_samples_nominal = Field(-1, "Nominal number of samples of the waveforms")
    num_samples_long = Field(-1, "Number of samples of the long waveforms")
    num_samples_removed_start = Field(-1, "Number of samples removed at the start of the waveforms")
    num_samples_removed_end = Field(-1, "Number of samples removed at the end of the waveforms")
    sampling_frequency = Field(-1, "Sampling frequency (MHz), DL0 only")
    data_model_version = Field("", "Version of the data model")
    calibration_service_id = Field(-1, "Id of the calibration service")
    calibration_algorithm_id = Field(-1, "Id of the calibration algorithm")


class SST1MObservationBlockContainer(ObservationBlockContainer):
    """
    ctapipe observation block of a SST-1M run (one raw data file), with the target
    of the TARGET field of the file header
    """

    default_prefix = ""

    target = Field("", "Target of the run, e.g. Crab, Transition or dark (TARGET field of the file)", max_length=64)
    wobble = Field(
        "NONE",
        "Wobble of the run, e.g. W1: UNDEF if the TARGET field has no wobble, NONE if it only has the target",
        max_length=16,
    )


class SST1MArrayEventContainer(ArrayEventContainer):
    """
    ctapipe array event with the SST-1M raw data (DigiCam specific R0 fields)
    and the SST-1M monitoring (statistics of the ADC samples of the pedestal events)
    """
    r0 = Field(default_factory=SST1MR0Container, description="Raw data of the SST-1M telescopes")
    monitoring = Field(default_factory=SST1MMonitoringContainer, description="Monitoring data")


class DL1_info(IsDescription):
    sst1mpipe_version    = StringCol(68)
    target               = StringCol(256)
    ra                   = Float64Col()
    dec                  = Float64Col()
    wobble               = StringCol(68)
    manual_coords        = BoolCol()
    calib_file           = StringCol(256)
    window_file          = StringCol(256)
    n_saturated          = Int64Col()
    n_pedestal           = Int64Col()
    n_survived_pedestals = Float64Col()
    raised_threshold_fraction = Float64Col()
    n_triggered_tel1     = Int64Col()
    n_triggered_tel2     = Int64Col()
    swat_event_ids_used  = BoolCol()


class DL2_info(IsDescription):
    sst1mpipe_version   = StringCol(68)
    RF_used             = StringCol(256)
