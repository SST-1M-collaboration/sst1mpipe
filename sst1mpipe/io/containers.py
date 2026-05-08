"""
Container structures for data that should be read or written to disk. The main
data container is DataContainer() and holds the containers of each data
processing level. The data processing levels start from R0 up to DL2, where R0
holds the cameras raw data and DL2 the air shower high-level parameters.
In general each major pipeline step is associated with a given data level.
Please keep in mind that the data level definition and the associated fields
might change rapidly as there is no final data level definition.
"""
# from aenum import IntFlag
from enum import Flag

#import pickle
#from gzip import open as gzip_open
#from os import remove
#from os.path import isfile
import numpy as np
from astropy import units as u
from ctapipe.containers import (
    ArrayEventContainer,
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
           'InstrumentContainer',
           'R0Container',
           'R0CameraContainer',
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


class InstrumentContainer(Container):
    """Storage of header info that does not change with event. This is a
    temporary hack until the Instrument module and database is fully
    implemented.  Eventually static information like this will not be
    part of the data stream, but be loaded and accessed from
    functions.
    """

    subarray = Field(None,
                     "SubarrayDescription from the instrument module")
    optical_foclen = Field(Map(np.ndarray), "map of tel_id to focal length")
    tel_pos = Field(Map(np.ndarray), "map of tel_id to telescope position")
    pixel_pos = Field(Map(np.ndarray), "map of tel_id to pixel positions")
    telescope_ids = Field([], "list of IDs of telescopes used in the run")
    num_pixels = Field(Map(int), "map of tel_id to number of pixels in camera")
    num_channels = Field(Map(int), "map of tel_id to number of channels")
    num_samples = Field(Map(int), "map of tel_id to number of samples")
    geom = Field(Map(None), 'map of tel_if to CameraGeometry')
    cam = Field(Map(None), 'map of tel_id to Camera')
    optics = Field(Map(None), 'map of tel_id to CameraOptics')
    cluster_matrix_7 = Field(Map(np.ndarray), 'map of tel_id of cluster 7 matrix')
    cluster_matrix_19 = Field(
        Map(np.ndarray),
        'map of tel_id of cluster 19 matrix'
    )
    patch_matrix = Field(Map(np.ndarray), 'map of tel_id of patch matrix')
    mirror_dish_area = Field(Map(float),
                             "map of tel_id to the area of the mirror dish",
                             unit=u.m ** 2)
    mirror_numtiles = Field(Map(int),
                            "map of tel_id to the number of \
                            tiles for the mirror")


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


class R0CameraContainer(Container):
    """
    Storage of raw data from a single telescope

    :param pixel_flags: a ndarray to flag the pixels
    :type pixel_flags: ndarray (n_pixels, ) (bool)
    :param adc_samples: a ndarray (n_pixels, n_samples) containing
                        the waveforms in each pixel
    :type adc_samples: ndarray (n_pixels, n_samples, ) (uint16)
    :param adc_sums: numpy array containing integrated ADC data
                     (n_channels, x n_pixels)
    :type adc_sums: ndarray (n_channels, x n_pixels)
    :param baseline: baseline holder for baseline computation using
                     clocked triggers
    :type baseline: ndarray (n_pixels, ) (float)
    :param digicam_baseline: baseline computed by DigiCam of pre-samples
                             (using 1024 samples)
    :type digicam_baseline: ndarray (n_pixels, ) (uint16)
    :param standard_deviation: baseline standard deviation holder for baseline
                               computed using clocked triggers
    :type standard deviation: ndarray (n_pixels, ) (float)
    :param dark_baseline: baseline holder for baseline computed in dark
                          condition (lid closed)
    :type dark_baseline: ndarray (n_pixels, ) (float)
    :param hv_off_baseline: baseline computed with sensors just bellow
                            breakdown voltage (or without bias voltage applied)
    :type hv_off_baseline: ndarray (n_pixels, ) (float)
    :param camera_event_id: unique event identification provided by DigiCam
    :type camera_event_id: (int)
    :param camera_event_number: event number within the first trigger of
                                operation
    :type camera_event_number: (int)
    :param local_camera_clock: time stamp from internal DigiCam clock (ns)
    :type local_camera_clock: (int)
    :param gps_time: time stamp provided by a precise external clock
                     (synchronized between hardware components)
    :type gps_time: (int)
    """
    pixel_flags = Field(np.ndarray, 'numpy array containing pixel flags')
    adc_samples = Field(np.ndarray,
                        "numpy array containing ADC samples"
                        "(n_channels x n_pixels, n_samples)")
    adc_sums = Field(np.ndarray, "numpy array containing integrated ADC data"
                              "(n_channels, x n_pixels)")
    baseline = Field(None, "number of time samples for telescope")
    digicam_baseline = Field(np.ndarray, 'Baseline computed by DigiCam')
    standard_deviation = Field(np.ndarray, "number of time samples for telescope")
    dark_baseline = Field(np.ndarray, 'dark baseline')
    hv_off_baseline = Field(np.ndarray, 'HV off baseline')
    camera_event_id = Field(int, 'Camera event number')
    camera_event_number = Field(int, "camera event number")
    local_camera_clock = Field(np.int64, "camera timestamp")
    gps_time = Field(np.int64, "gps timestamp")
    white_rabbit_time = Field(float, "precise white rabbit based timestamp")
    _camera_event_type = Field(CameraEventType, "camera event type")

    @property
    def camera_event_type(self):
        return self._camera_event_type

    @camera_event_type.setter
    def camera_event_type(self, value):
        self._camera_event_type = CameraEventType(value)

    array_event_type = Field(int, "array event type")
    trigger_input_traces = Field(np.ndarray, "trigger patch trace (n_patches)")
    trigger_input_offline = Field(np.ndarray, "trigger patch trace (n_patches)")
    trigger_output_patch7 = Field(np.ndarray, "trigger 7 patch cluster trace \
                                  (n_clusters)")
    trigger_output_patch19 = Field(np.ndarray, "trigger 19 patch cluster trace \
                                   (n_clusters)")
    trigger_input_7 = Field(np.ndarray, 'trigger input CLUSTER7')
    trigger_input_19 = Field(np.ndarray, 'trigger input CLUSTER19')
    num_samples = Field(int, "number of time samples for telescope")


class R0Container(Container):
    """
    Storage of a Merged Raw Data Event
    """

    run_id = Field(-1, "run id number")
    event_id = Field(-1, "event id number")
    tels_with_data = Field([], "list of telescopes with data")
    tel = Field(Map(R0CameraContainer), "map of tel_id to R0CameraContainer")


class R1CameraContainer(Container):
    """
    Storage of r1 calibrated data from a single telescope
    """

    adc_samples = Field(np.ndarray, "baseline subtracted ADCs, (n_pixels, \
                        n_samples)")
    nsb = Field(np.ndarray, "nsb rate in GHz")
    pde = Field(np.ndarray, "Photo Detection Efficiency at given NSB")
    gain_drop = Field(np.ndarray, "gain drop")


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


class SST1MContainer(Container):
    r0 = Field(R0Container(), "Raw Data")
    r1 = Field(R1Container(), "Raw Common Data")

    inst = Field(InstrumentContainer(), "Instrumental information")
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
    sst1m = Field(SST1MContainer(), "SST1M specific information")



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
