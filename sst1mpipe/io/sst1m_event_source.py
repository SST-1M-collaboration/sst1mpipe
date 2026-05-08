
import warnings

import numpy as np
from astropy import units as u
from astropy.time import Time
from ctapipe.containers import (
    ObservationBlockContainer,
    PointingMode,
    SchedulingBlockContainer,
)
from ctapipe.core.traits import Bool, Float
from ctapipe.instrument import (
    SubarrayDescription,
)
from ctapipe.instrument.subarray import EarthLocation
from ctapipe.io import (
    EventSource,
)
from ctapipe.io.datalevels import DataLevel
from protozfits import File, MultiZFitsFiles

from sst1mpipe.constants import (
    PATCH_ID_INPUT_SORT_IDS,
    PATCH_ID_OUTPUT_SORT_IDS,
    REFERENCE_LOCATION,
)
from sst1mpipe.instrument import camera
from sst1mpipe.io.containers import (
    SST1MArrayEventContainer,
)

# from tqdm import tqdm

# from sst1mpipe.io.zfits import (
#     _prepare_trigger_input,
#     _prepare_trigger_output
# )




class SST1MEventSource(EventSource):
    """
    https://github.com/cta-observatory/ctapipe_io_lst/blob/0f8b8cd39403f51dc8b1b0e1eb5a6045ea5deb15/src/ctapipe_io_lst/__init__.py#L13
    EventSource for SST1M R0 data.

    Reimplementation of CTAPIPE_IO_LST for SST1M which also uses fits input files which are not readable by a default ctapipe
    """

    reference_position_lon = Float(
        default_value = REFERENCE_LOCATION.lon.deg,
        help = (
            "Longitude of the reference location for telescope GroundFrame coordinates."
        )
    ).tag(config = True)

    reference_position_lat = Float(
        default_value = REFERENCE_LOCATION.lat.deg,
        help = (
            "Latitude of the reference location for telescope GroundFrame coordinates."
        )
    ).tag(config = True)

    reference_position_height = Float(
        default_value = REFERENCE_LOCATION.height.to_value(u.m),
        help = (
            "Height of the reference location for telescope GroundFrame coordinates."
            " Default is current MC obslevel."
        )
    ).tag(config = True)

    pointing_information = Bool(
        default_value = True,
        help = (
            'Fill pointing information.'
        ),
    ).tag(config = True)

    def __init__(self,
                 filelist=None,
                 camera=camera.DigiCam,
                 max_events=None,
                 event_id = None,
                 allowed_tels = None,
                 disable_bar = False,
                 **kwargs
        ):
        # LST/CTA uses differenct filename naming convention, how to work with the SST1M file naming convention?
        # for nowadays EventSource obtains only the first file,
        # but SST1MEventSource counts with all input files, implemented via MultiZFitsFiles

        super().__init__(input_url=filelist[0], **kwargs)

        self.filelist = filelist
        self.run_id = 0
        self.tel_id = 0
        self.camera = camera
        self.max_events = max_events
        self.event_id = event_id
        self.allowed_tels = allowed_tels
        self.disable_bar = disable_bar

        # LST reads camera_config from input files, is it needed such functionality for SST1M?
        self.camera_config = None
        self.run_start = Time(self.camera_config.date, format='unix') if self.camera_config is not None else None

        reference_location = EarthLocation(
            lon = self.reference_position_lon,
            lat = self.reference_position_lat,
            height = self.reference_position_height,
        )
        self._subarray = self.create_subarray(self.tel_id, reference_location)


        # self.pointing_source = PointingSource(subarray=self.subarray, parent=self)

        target_info = {}
        pointing_mode = PointingMode.UNKNOWN
        # if self.pointing_information:
        #     target = self.pointing_source.get_target(tel_id=self.tel_id, time=self.run_start)
        #     if target is not None:
        #         target_info["subarray_pointing_lon"] = target["ra"]
        #         target_info["subarray_pointing_lat"] = target["dec"]
        #         target_info["subarray_pointing_frame"] = CoordinateFrameType.ICRS
        #         pointing_mode = PointingMode.TRACK

        self._scheduling_blocks = {
            self.run_id: SchedulingBlockContainer(
                sb_id=np.uint64(self.run_id),
                producer_id=f"SST1M-{self.tel_id}",
                pointing_mode=pointing_mode,
            )
        }

        self._observation_blocks = {
            self.run_id: ObservationBlockContainer(
                obs_id=np.uint64(self.run_id),
                sb_id=np.uint64(self.run_id),
                producer_id=f"SST1M-{self.tel_id}",
                actual_start_time=self.run_start,
                **target_info
            )
        }

        self._swat_event_ids_available = self.check_swat_event_ids_available(filelist)

    @property
    def subarray(self):
        return self._subarray

    @property
    def is_simulation(self):
        return False

    # @property
    # def obs_ids(self):
    #     # currently no obs id is available from the input files
    #     return list(self.observation_blocks)

    @property
    def observation_blocks(self):
        return self._observation_blocks

    @property
    def scheduling_blocks(self):
        return self._scheduling_blocks

    @property
    def datalevels(self):
        # if self.r0_r1_calibrator.calibration_path is not None:
        #     return (DataLevel.R0, DataLevel.R1)
        return (DataLevel.R0, )

    @property
    def swat_event_ids_available(self):
        return self._swat_event_ids_available

    @staticmethod
    def check_swat_event_ids_available(filelist):
        """
        Determine if the files contain SWAT-generated arrayEvtNum IDs
        Returns
        -------
        True  if the files contain array-level IDs
        False if the files do not contain array-level IDs
        """
        with File(filelist[0]) as f:
            id0 = f.Events[0].arrayEvtNum
            id1 = f.Events[1].arrayEvtNum
        # If SWAT arrayEvtNum was not written, the value is always 0
        # Otherwise, we can expect 0 at most once, if SWAT was just restarted
        if id0 == id1 == 0:
            return False
        else:
            return True

    @staticmethod
    def create_subarray(tel_id=1, reference_location=None):
        """
        Obtain the subarray from the EventSource
        Returns
        -------
        ctapipe.instrument.SubarrayDescription
        """

        subarray = SubarrayDescription(
            name=f"SST1M-{tel_id} subarray",
            tel_positions={},
            tel_descriptions={},
            reference_location=reference_location,
        )

        return subarray

    def _generator(self):
        """
        """
        yield from self.get_array_event(self.filelist)

    def get_array_event(self, input_path : str):
        """
        """
        print(f"input_path : {input_path}")
        loaded_telescopes = []
        array_event = SST1MArrayEventContainer()
        with MultiZFitsFiles(input_path) as events:
            array_event.r0.meta = dict(is_simulation=False)
            for event_counter, event in enumerate(events):
                # print(f" **** event: {event}")
                if self.max_events is not None and event_counter > self.max_events:
                    break
                array_event.count = event_counter
                if self._swat_event_ids_available:
                    array_event.sst1m.r0.event_id = event.arrayEvtNum
                else:
                    array_event.sst1m.r0.event_id = event.eventNumber
                array_event.sst1m.r0.tels_with_data = [event.telescopeID, ]
                _sort_ids = None
                for tel_id in array_event.sst1m.r0.tels_with_data:
                    pixel_ids = event.hiGain.waveforms.pixelsIndices
                    n_pixels = len(pixel_ids)
                    if _sort_ids is None:
                        _sort_ids = np.argsort(pixel_ids)
                    samples = event.hiGain.waveforms.samples.reshape(n_pixels, -1)

                    try:
                        unsorted_baseline = event.hiGain.waveforms.baselines
                    except AttributeError:
                        warnings.warn(
                            "Could not read `hiGain.waveforms.baselines`"
                            f"for event:{event_counter} (eventNumber {event.eventNumber})\n"
                            f"of file:{self.input_url}\n"
                        )
                        return np.ones(n_pixels) * np.nan

                    if tel_id not in loaded_telescopes:
                        array_event.sst1m.inst.num_channels[tel_id] = event.num_gains
                        array_event.sst1m.inst.geom[tel_id] = self.camera.geometry
                        array_event.sst1m.inst.cluster_matrix_7[tel_id] = \
                            self.camera.cluster_7_matrix
                        array_event.sst1m.inst.cluster_matrix_19[tel_id] = \
                            self.camera.cluster_19_matrix
                        array_event.sst1m.inst.patch_matrix[tel_id] = self.camera.patch_matrix
                        array_event.sst1m.inst.num_pixels[tel_id] = samples.shape[0]
                        array_event.sst1m.inst.num_samples[tel_id] = samples.shape[1]
                        loaded_telescopes.append(tel_id)

                    cta_r0 = array_event.r0.tel[tel_id]
                    cta_r0.waveform = samples[_sort_ids].reshape(1, n_pixels, -1)

                    r0 = array_event.sst1m.r0.tel[tel_id]
                    r0.camera_event_number = event.eventNumber
                    r0.pixel_flags = event.pixels_flags[_sort_ids]
                    r0.local_camera_clock = (
                        np.int64(event.local_time_sec * 1E9) +
                        np.int64(event.local_time_nanosec)
                    )
                    r0.gps_time = (
                        np.int64(event.trig.timeSec * 1E9) +
                        np.int64(event.trig.timeNanoSec)
                    )
                    r0.camera_event_type = event.event_type
                    r0.array_event_type = event.eventType
                    r0.adc_samples = samples[_sort_ids]

                    if len(event.trigger_input_traces) > 0:
                        r0.trigger_input_traces = self._prepare_trigger_input(
                            event.trigger_input_traces
                        )
                    else:
                        warnings.warn(
                            'trigger_input_traces does not exist: --> nan')
                        r0.trigger_input_traces = np.zeros(
                            (432, array_event.sst1m.inst.num_samples[tel_id])) * np.nan

                    if len(event.trigger_output_patch7) > 0:
                        r0.trigger_output_patch7 = self._prepare_trigger_output(
                            event.trigger_output_patch7)
                    else:
                        warnings.warn(
                            'trigger_output_patch7 does not exist: --> nan')
                        r0.trigger_output_patch7 = np.zeros(
                            (432, array_event.sst1m.inst.num_samples[tel_id])) * np.nan

                    if len(event.trigger_output_patch19) > 0:
                        r0.trigger_output_patch19 = self._prepare_trigger_output(
                            event.trigger_output_patch19)
                    else:
                        warnings.warn(
                            'trigger_output_patch19 does not exist: --> nan')
                        r0.trigger_output_patch19 = np.zeros(
                            (432, array_event.sst1m.inst.num_samples[tel_id])) * np.nan

                    r0.digicam_baseline = unsorted_baseline[_sort_ids] / 16
                yield array_event

    def _prepare_trigger_input(self, _a):
        A, B = 3, 192
        cut = 144
        _a = _a.reshape(-1, A)
        _a = _a.reshape(-1, A, B)
        _a = _a[..., :cut]
        _a = _a.reshape(_a.shape[0], -1)
        _a = _a.T
        _a = _a[PATCH_ID_INPUT_SORT_IDS]
        return _a


    def _prepare_trigger_output(self, _a):
        A, B, C = 3, 18, 8

        _a = np.unpackbits(_a.reshape(-1, A, B, 1), axis=-1)
        _a = _a[..., ::-1]
        _a = _a.reshape(-1, A * B * C).T
        return _a[PATCH_ID_OUTPUT_SORT_IDS]

    @staticmethod
    def is_compatible(file_path):
        pass
        # from astropy.io import fits

        # try:
        #     with fits.open(file_path) as hdul:
        #         if "Events" not in hdul:
        #             return False

        #         header = hdul["Events"].header
        #         ttypes = {
        #             value for key, value in header.items()
        #             if 'TTYPE' in key
        #         }
        # except OSError:
        #     return False


        # is_protobuf_zfits_file = (
        #     (header['XTENSION'] == 'BINTABLE')
        #     and (header['ZTABLE'] is True)
        #     and (header['ORIGIN'] == 'CTA')
        #     and (header['PBFHEAD'] == 'R1.CameraEvent')
        # )

        # print(header["XTENSION"], header["ZTABLE"], header["ORIGIN"], header["PBFHEAD"])
        # return True
        # # is_lst_file = 'lstcam_counters' in ttypes
        # # return is_protobuf_zfits_file & is_lst_file
