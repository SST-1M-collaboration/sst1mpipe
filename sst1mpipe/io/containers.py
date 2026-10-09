"""
Containers of the SST-1M data. The events are ctapipe `~ctapipe.containers.ArrayEventContainer`
(`SST1MArrayEventContainer`) with the DigiCam specific raw data (R0) and the statistics of
the ADC samples of the pedestal events (monitoring). The other data levels follow ctapipe.
The PyTables descriptions are those of the tables written by sst1mpipe in the DL1/DL2 files.
"""
from enum import Flag
from functools import partial

import astropy.units as u
import numpy as np
from ctapipe.containers import (
    ArrayEventContainer,
    MonitoringCameraContainer,
    MonitoringContainer,
    ObservationBlockContainer,
    NAN_TIME,
    PedestalContainer,
    R0CameraContainer,
    R0Container,
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
    "SST1MR0CameraContainer",
    "SST1MR0Container",
    "R0PedestalContainer",
    "TelescopePointingMonitoringContainer",
    "SST1MMonitoringCameraContainer",
    "SST1MMonitoringContainer",
    "SST1MArrayEventContainer",
    "SST1MObservationBlockContainer",
    "DigicamConfigContainer",
    "StreamCameraConfigContainer",
    "DL1_info",
    "DL2_info",
]


class CameraEventType(Flag):
    # from https://github.com/cta-sst-1m/digicampipe/issues/244
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
    pedestal = Field(None, "baseline computed by DigiCam from 1024 pre-samples, in ADC (n_pixels)")
    camera_event_number = Field(None, "event number within the first trigger of operation")
    event_time = Field(None, "timestamp from the internal DigiCam clock (ns, TAI)")
    gps_time = Field(None, "timestamp from a precise external clock (ns)")
    white_rabbit_time = Field(None, "precise White Rabbit based timestamp")
    _event_type = Field(None, "camera event type")
    trigger_input_traces = Field(None, "trigger patch traces (n_patches, n_samples)")
    trigger_output_patch7 = Field(None, "trigger 7 patch cluster traces (n_clusters, n_samples)")
    trigger_output_patch19 = Field(None, "trigger 19 patch cluster traces (n_clusters, n_samples)")
    trigger_output_muon = Field(None, "trigger muon cluster traces (n_clusters, n_samples)")

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


class R0PedestalContainer(PedestalContainer):
    """
    Statistics of the ADC samples (in ADC) of the pedestal events.
    The statistics of the calibrated images (in p.e.) are stored in the ctapipe
    `~ctapipe.containers.PedestalContainer` of `~ctapipe.containers.MonitoringCameraContainer`.
    """

    default_prefix = "pedestal"


class TelescopePointingMonitoringContainer(Container):
    """
    Pointing of the telescope (alt/az) at a time, written in the table
    ``/dl0/monitoring/telescope/pointing/tel_XXX`` read by the ctapipe
    ``PointingInterpolator`` (``HDF5EventSource``, ``TableLoader``)
    """

    default_prefix = ""

    time = Field(NAN_TIME, "Time of the pointing")
    azimuth = Field(np.nan * u.rad, "Azimuth of the pointing", unit=u.rad)
    altitude = Field(np.nan * u.rad, "Altitude of the pointing", unit=u.rad)


class SST1MMonitoringCameraContainer(MonitoringCameraContainer):
    """
    ctapipe camera monitoring with the R0 level monitoring and the pointing of SST-1M
    """

    r0 = Field(
        default_factory=R0PedestalContainer,
        description="Statistics of the ADC samples of the pedestal events",
    )
    pointing = Field(
        default_factory=TelescopePointingMonitoringContainer,
        description="Pointing of the telescope at the time of the event",
    )


class SST1MMonitoringContainer(MonitoringContainer):
    tel = Field(
        default_factory=partial(Map, SST1MMonitoringCameraContainer),
        description="map of tel_id to SST1MMonitoringCameraContainer",
    )


class DigicamConfigContainer(Container):
    """
    Configuration of the DigiCam boards of the camera (``DigicamConfig`` table of the
    raw data file). The arrays have one entry per board slot (n_boards, 39), 0 for
    the empty slots.
    """

    default_prefix = "digicam"

    protocol_vers = Field(None, "Protocol version of each board (n_boards)")
    sn = Field(None, "Serial number of each board (n_boards)")
    hv = Field(None, "High voltage status of each board (n_boards)")
    gateware_rev = Field(None, "Gateware revision of each board (n_boards)")
    gateware_vers = Field(None, "Gateware version of each board (n_boards)")
    gateware_code = Field(None, "Gateware code of each board (n_boards)")
    gateware_card_type = Field(None, "Card type of the gateware of each board (n_boards)")
    firmware_rev = Field(None, "Firmware revision of each board (n_boards)")
    firmware_vers = Field(None, "Firmware version of each board (n_boards)")
    firmware_code = Field(None, "Firmware code of each board (n_boards)")
    firmware_card_type = Field(None, "Card type of the firmware of each board (n_boards)")
    operation_id = Field(-1, "Id of the operation of DigiCam")
    operation_data = Field(-1, "Data of the operation of DigiCam")
    digicam_time_sec = Field(-1, "DigiCam time of the configuration, seconds")
    digicam_time_nanosec = Field(-1, "DigiCam time of the configuration, nanoseconds")


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
    and the R0 level monitoring (statistics of the ADC samples of the pedestal events)
    """
    r0 = Field(default_factory=SST1MR0Container, description="Raw data of the SST-1M telescopes")
    mon = Field(default_factory=SST1MMonitoringContainer, description="container for monitoring data (MON)")


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
    n_triggered_tel1     = Int64Col()
    n_triggered_tel2     = Int64Col()
    swat_event_ids_used  = BoolCol()


class DL2_info(IsDescription):
    sst1mpipe_version   = StringCol(68)
    RF_used             = StringCol(256)
