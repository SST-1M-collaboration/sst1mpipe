"""
Container structures for data that should be read or written to disk. The main
data container is DataContainer() and holds the containers of each data
processing level. The data processing levels start from R0 up to DL2, where R0
holds the cameras raw data and DL2 the air shower high-level parameters.
In general each major pipeline step is associated with a given data level.
Please keep in mind that the data level definition and the associated fields
might change rapidly as there is no final data level definition.
"""
from enum import Flag
from functools import partial

import numpy as np
from ctapipe.containers import (
    ArrayEventContainer,
    MonitoringCameraContainer,
    MonitoringContainer,
    PedestalContainer,
    R0CameraContainer,
    R0Container,
    TriggerContainer,
)
from ctapipe.core import Container, Field, Map

# from ctapipe.serializer import Serializer
# from ctapipe.containers import MCEventContainer, ReconstructedContainer, \
# MCHeaderContainer, CentralTriggerContainer
from tables import (
    BoolCol,
    Float64Col,
    Int64Col,
    IsDescription,
    StringCol,
)

__all__ = ['CameraEventType',
           'SST1MR0Container',
           'SST1MR0CameraContainer',
           'R1Container',
           'R1CameraContainer',
        #    'DL0Container',
        #    'DL0CameraContainer',
        #    'DL1Container',
        #    'DL1CameraContainer',
        #    'MCEventContainer',
        #    'DataContainer',
            "SST1MArrayEventContainer",
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


# class DL1CameraContainer(Container):
#     """Storage of output of camera calibrationm e.g the final calibrated
#     image in intensity units and other per-event calculated
#     calibration information.
#     """

#     pe_samples = Field(np.ndarray, "numpy array containing data volume reduced \
#                        p.e. samples (n_channels x n_pixels)")
#     cleaning_mask = Field(np.ndarray, "mask for clean pixels")
#     time_bin = Field(np.ndarray, "numpy array containing the bin of maximum \
#                     (n_pixels)")
#     pe_samples_trace = Field(np.ndarray, "numpy array containing data volume \
#                              reduced p.e. samples (n_channels x n_pixels, \
#                              n_samples)")
#     on_border = Field(bool, "Boolean telling if the shower touches the camera \
#                       border or not")
#     time_spread = Field(float, 'Time elongation of the shower')


# class DL1Container(Container):
#     """ DL1 Calibrated Camera Images and associated data"""
#     tel = Field(Map(DL1CameraContainer), "map of tel_id to DL1CameraContainer")


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
    adc_samples = Field(
        None,
        "ADC samples (n_pixels, n_samples). Same as waveform[0] when read, but a separate"
        " array: the waveform of the bad pixels is set to 0 by the calibration, not adc_samples",
    )
    adc_sums = Field(None, "numpy array containing integrated ADC data (n_channels, n_pixels)")
    baseline = Field(None, "baseline computed using clocked triggers (n_pixels)")
    digicam_baseline = Field(None, "baseline computed by DigiCam from 1024 pre-samples (n_pixels)")
    standard_deviation = Field(None, "baseline standard deviation computed using clocked triggers (n_pixels)")
    dark_baseline = Field(None, "baseline computed in dark condition, lid closed (n_pixels)")
    hv_off_baseline = Field(None, "baseline computed without bias voltage (n_pixels)")
    camera_event_id = Field(None, "unique event identification provided by DigiCam")
    camera_event_number = Field(None, "event number within the first trigger of operation")
    local_camera_clock = Field(None, "timestamp from the internal DigiCam clock (ns, TAI)")
    gps_time = Field(None, "timestamp from a precise external clock (ns)")
    white_rabbit_time = Field(None, "precise White Rabbit based timestamp")
    _camera_event_type = Field(None, "camera event type")
    array_event_type = Field(None, "array event type")
    trigger_input_traces = Field(None, "trigger patch traces (n_patches, n_samples)")
    trigger_input_offline = Field(None, "trigger patch traces computed offline (n_patches, n_samples)")
    trigger_output_patch7 = Field(None, "trigger 7 patch cluster traces (n_clusters, n_samples)")
    trigger_output_patch19 = Field(None, "trigger 19 patch cluster traces (n_clusters, n_samples)")
    trigger_input_7 = Field(None, "trigger input CLUSTER7")
    trigger_input_19 = Field(None, "trigger input CLUSTER19")
    num_samples = Field(None, "number of time samples")

    @property
    def camera_event_type(self):
        return self._camera_event_type

    @camera_event_type.setter
    def camera_event_type(self, value):
        self._camera_event_type = CameraEventType(value)


class SST1MR0Container(R0Container):
    """
    Raw data of the SST-1M telescopes
    """

    tel = Field(
        default_factory=partial(Map, SST1MR0CameraContainer),
        description="map of tel_id to SST1MR0CameraContainer",
    )


class R1CameraContainer(Container):
    """
    Storage of r1 calibrated data from a single telescope
    """

    adc_samples = Field(np.ndarray, "baseline subtracted ADCs, (n_pixels, \
                        n_samples)")
    nsb = Field(np.ndarray, "nsb rate in GHz")
    pde = Field(np.ndarray, "Photo Detection Efficiency at given NSB")
    gain_drop = Field(np.ndarray, "gain drop")
    saturated = Field(False, "True if the charge of saturated pixels was corrected")


class R1Container(Container):
    """
    Storage of a r1 calibrated Data Event
    """

    run_id = Field(-1, "run id number")
    event_id = Field(-1, "event id number")
    tels_with_data = Field([], "list of telescopes with data")
    tel = Field(Map(R1CameraContainer), "map of tel_id to R1CameraContainer")


# class DL0CameraContainer(Container):
#     """
#     Storage of data volume reduced dl0 data from a single telescope
#     """


# class DL0Container(Container):
#     """
#     Storage of a data volume reduced Event
#     """

#     run_id = Field(-1, "run id number")
#     event_id = Field(-1, "event id number")
#     tels_with_data = Field([], "list of telescopes with data")
#     tel = Field(Map(DL0CameraContainer), "map of tel_id to DL0CameraContainer")


class R0PedestalContainer(PedestalContainer):
    """
    Statistics of the ADC samples (in ADC) of the pedestal events.
    The statistics of the calibrated images (in p.e.) are stored in the ctapipe
    `~ctapipe.containers.PedestalContainer` of `~ctapipe.containers.MonitoringCameraContainer`.
    """

    default_prefix = "pedestal"


class SST1MMonitoringCameraContainer(MonitoringCameraContainer):
    """
    ctapipe camera monitoring with the R0 level monitoring of SST-1M
    """

    r0 = Field(
        default_factory=R0PedestalContainer,
        description="Statistics of the ADC samples of the pedestal events",
    )


class SST1MMonitoringContainer(MonitoringContainer):
    tel = Field(
        default_factory=partial(Map, SST1MMonitoringCameraContainer),
        description="map of tel_id to SST1MMonitoringCameraContainer",
    )


class SST1MContainer(Container):
    r1 = Field(R1Container(), "SST-1M specific information of the calibration")
    slow_data = Field(None, "Slow Data Information")
    trig = Field(TriggerContainer(), "central trigger information")
    count = Field(0, "number of events processed")
    # dl0 = Field(DL0Container(), "DL0 Data Volume Reduced Data")
    # dl1 = Field(DL1Container(), "DL1 Calibrated image")
    # dl2 = Field(ReconstructedContainer(), "Reconstructed Shower Information")
    # mc = Field(MCEventContainer(), "Monte-Carlo data")
    # mcheader = Field(MCHeaderContainer(), "Monte-Carlo run header data")

class SST1MArrayEventContainer(ArrayEventContainer):
    """
    Data container including SST1M and monitoring information
    """
    r0 = Field(default_factory=SST1MR0Container, description="Raw data of the SST-1M telescopes")
    sst1m = Field(SST1MContainer(), "SST1M specific information")
    mon = Field(default_factory=SST1MMonitoringContainer, description="container for monitoring data (MON)")



# class DataContainer(Container):
#     """ Top-level container for all event information.
#     Each field is representing a specific data processing level from (R0 to
#     DL2) Please keep in mind that the data level definition and the associated
#     fields might change rapidly as there is not a final data format. The data
#     levels R0, R1, DL1, contains sub-containers for each telescope.
#     After DL2 the data is not processed at the telescope level.
#     """
#     r0 = Field(R0Container(), "Raw Data")
#     r1 = Field(R1Container(), "Raw Common Data")
#     # dl0 = Field(DL0Container(), "DL0 Data Volume Reduced Data")
#     # dl1 = Field(DL1Container(), "DL1 Calibrated image")
#     # dl2 = Field(ReconstructedContainer(), "Reconstructed Shower Information")
#     # mc = Field(MCEventContainer(), "Monte-Carlo data")
#     # mcheader = Field(MCHeaderContainer(), "Monte-Carlo run header data")
#     # inst = Field(InstrumentContainer(), "Instrumental information")
#     slow_data = Field(None, "Slow Data Information")
#     # trig = Field(CentralTriggerContainer(), "central trigger information")
#     count = Field(0, "number of events processed")


# def load_from_pickle_gz(file):
#     file = gzip_open(file, "rb")
#     while True:
#         try:
#             yield pickle.load(file)
#         except (EOFError, pickle.UnpicklingError):
#             return


# def save_to_pickle_gz(event_stream, file, overwrite=False, max_events=None):
#     if isfile(file):
#         if overwrite:
#             print('remove old', file, 'file')
#             remove(file)
#         else:
#             print(file, 'exist, exiting...')
#             return
#     writer = Serializer(filename=file, format='pickle', mode='w')
#     counter_events = 0
#     for event in event_stream:
#         writer.add_container(event)
#         counter_events += 1

#         if max_events is not None and counter_events >= max_events:
#             break

#     writer.close()


# class CalibrationEventContainer(Container):
#     """
#     description test
#     """
#     # Raw

#     adc_samples = Field(np.ndarray, 'the raw data')
#     digicam_baseline = Field(np.ndarray, 'the baseline computed by the camera')
#     local_time = Field(np.ndarray, 'timestamps')
#     gps_time = Field(np.ndarray, 'time')

#     # Processed

#     dark_baseline = Field(np.ndarray, 'the baseline computed in dark')
#     baseline_shift = Field(np.ndarray, 'the baseline shift')
#     nsb_rate = Field(np.ndarray, 'Night sky background rate')
#     gain_drop = Field(np.ndarray, 'Gain drop')
#     baseline = Field(np.ndarray, 'the reconstructed baseline')
#     baseline_std = Field(np.ndarray, 'Baseline std')
#     pulse_mask = Field(np.ndarray, 'mask of adc_samples. True if the adc sample'
#                                 'contains a pulse  else False')
#     reconstructed_amplitude = Field(np.ndarray, 'array of the same shape as '
#                                              'adc_samples giving the'
#                                              ' reconstructed pulse amplitude'
#                                              ' for each adc sample')
#     reconstructed_charge = Field(np.ndarray, 'array of the same shape as '
#                                           'adc_samples giving the '
#                                           'reconstructed charge for each adc '
#                                           'sample')
#     reconstructed_number_of_pe = Field(np.ndarray, 'estimated number of photon '
#                                                 'electrons for each adc sample'
#                                        )
#     sample_pe = Field(
#         ndarray,
#         'array of the same shape as adc_samples giving the estimated fraction '
#         'of photon electrons for each adc sample'
#     )
#     reconstructed_time = Field(np.ndarray, 'reconstructed time '
#                                         'for each adc sample')
#     cleaning_mask = Field(np.ndarray, 'cleaning mask, pixel bool array')
#     shower = Field(bool, 'is the event considered as a shower')
#     border = Field(bool, 'is the event after cleaning touchin the camera '
#                          'borders')
#     burst = Field(bool, 'is the event during a burst')
#     saturated = Field(bool, 'is any pixel signal saturated')

#     def plot(self, pixel_id):
#         plt.figure()
#         plt.title('pixel : {}'.format(pixel_id))
#         plt.plot(self.adc_samples[pixel_id], label='raw')
#         plt.plot(self.pulse_mask[pixel_id], label='peak position')
#         plt.plot(self.reconstructed_charge[pixel_id], label='charge',
#                  linestyle='None', marker='o')
#         plt.plot(self.reconstructed_amplitude[pixel_id], label='amplitude',
#                  linestyle='None', marker='o')
#         plt.legend()


# class CalibrationContainerMeta(Container):
#     time = Field(float, 'time of the event')
#     event_id = Field(int, 'event id')
#     type = Field(int, 'event type')


# class CalibrationContainer(Container):
#     """
#     This Container() is used for the camera calibration pipeline.
#     It is meant to save each step of the calibration pipeline
#     """

#     config = Field(list, 'List of the input parameters'
#                          ' of the calibration analysis')  # Should use dict?
#     pixel_id = Field(np.ndarray, 'pixel ids')
#     data = CalibrationEventContainer()
#     event_id = Field(int, 'event_id')
#     event_type = Field(CameraEventType, 'Event type')
#     hillas = Field(HillasParametersContainer, 'Hillas parameters')
#     info = CalibrationContainerMeta()
#     slow_data = Field(None, "Slow Data Information")
#     mc = Field(MCEventContainer(), "Monte-Carlo data")


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
