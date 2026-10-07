
import logging
import os
import re
import warnings
from itertools import islice

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, SkyCoord
import astropy.io.ascii as aio
from astropy.io import fits
from astropy.time import Time
from ctapipe.containers import (
    CoordinateFrameType,
    EventType,
    ObservationBlockContainer,
    PointingMode,
    SchedulingBlockContainer,
)
from ctapipe.core import Provenance
from ctapipe.core.traits import Bool, Dict, Float, List, UseEnum
from ctapipe.instrument import FocalLengthKind
from ctapipe.io import (
    EventSource,
)
from ctapipe.io.datalevels import DataLevel
from protozfits import File

from sst1mpipe.constants import (
    PATCH_ID_INPUT_SORT_IDS,
    PATCH_ID_OUTPUT_SORT_IDS,
    REFERENCE_LOCATION,
    SUBARRAY_DESCRIPTION
)
from sst1mpipe.resources import PIXEL_MAPPING_FILE
from sst1mpipe.io.containers import (
    CameraEventType,
    SST1MArrayEventContainer,
)

logger = logging.getLogger(__name__)

# Number of events read at the beginning of each file to look for SWAT event ids
N_EVENTS_SWAT_ID_CHECK = 10


def parse_target_field(field):
    """
    Parse the TARGET field of the header of the ``Events`` table, expected as
    ``target[_,]wobble[_,]ra[_,]dec`` or ``target[_,]ra[_,]dec`` with ra, dec in deg,
    e.g. ``Crab_W1_83.63_22.01``. Files without pointing (e.g. ``dark`` or
    ``transition``) only contain the target.

    Parameters
    ----------
    field: str or None

    Returns
    -------
    target, wobble, ra, dec:
        None for the missing entries. wobble is None if the field has
        no delimiter and ``UNDEF`` if it has no ``W<n>`` entry.
    """
    if field is None:
        return None, None, None, None

    if field.count('_') > 1:
        delimiter = '_'
    elif field.count(',') > 1:
        delimiter = ','
    else:
        return field, None, None, None

    entries = field.split(delimiter)
    target = entries[0]
    match = re.search(r'W\d+', field)
    wobble = match.group(0) if match else 'UNDEF'

    if len(entries) not in (3, 4):
        return target, wobble, None, None
    try:
        ra, dec = float(entries[-2]), float(entries[-1])
    except ValueError:
        return target, wobble, None, None
    return target, wobble, ra, dec


def parse_file_name(file_name):
    """
    Date and run number of a SST-1M raw data file name ``SST1M<tel>_<date>_<run>.fits.fz``,
    e.g. ``SST1M1_20260121_0001.fits.fz`` -> ("20260121", "0001"). None if it does not match.
    """
    match = re.match(r'SST1M\d*_(\d+)_(\d+)', os.path.basename(str(file_name)))
    return match.groups() if match else None


def camera_clock_to_time(local_camera_clock):
    """
    Convert the camera clock (ns, TAI scale) to an astropy Time with ns precision, see
    https://github.com/cta-observatory/ctapipe_io_nectarcam/issues/24
    """
    localtime = np.uint64(local_camera_clock)
    S_TO_NS = np.uint64(1e9)
    full_seconds = localtime // S_TO_NS
    fractional_seconds = (localtime % S_TO_NS) / S_TO_NS
    return Time(full_seconds, fractional_seconds, format='unix_tai')


def file_has_swat_event_ids(path, n_events=N_EVENTS_SWAT_ID_CHECK):
    """
    True if any of the first ``n_events`` events of the file has a non zero
    array event id (``arrayEvtNum``) written by SWAT. False for an empty file.
    """
    with File(str(path)) as f:
        return any(event.arrayEvtNum != 0 for event in islice(f.Events, n_events))


def parse_swapped_modules(swapped_modules, pixel_mapping_file=PIXEL_MAPPING_FILE):
    """
    Pixels of the wrongly connected modules of the cameras, which must be swapped.

    Parameters
    ----------
    swapped_modules: list of dict
        Entries with the telescope (``tel_id``), the period (``start``, ``stop``,
        UTC, ISO format) and the two modules (``modules``)
    pixel_mapping_file: path
        DigiCam pixel mapping, with the module of each pixel (``pixel_sw_id``)

    Returns
    -------
    dict:
        tel_id -> list of (start, stop, pixels_1, pixels_2): period (unix TAI, s) and
        pixel ids of the two modules
    """
    if len(swapped_modules) == 0:
        return {}
    mapping = aio.read(pixel_mapping_file)

    def module_pixels(module):
        return np.sort(np.asarray(mapping[mapping["module"] == module]["pixel_sw_id"]))

    pixel_swaps = {}
    for entry in swapped_modules:
        start = Time(entry["start"], format="isot", scale="utc").unix_tai
        stop = Time(entry["stop"], format="isot", scale="utc").unix_tai
        module_1, module_2 = entry["modules"]
        pixel_swaps.setdefault(entry["tel_id"], []).append(
            (start, stop, module_pixels(module_1), module_pixels(module_2))
        )
    return pixel_swaps


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

    pointing_ra = Float(
        default_value=None,
        allow_none=True,
        help="Pointing right ascension in deg. Overrides the TARGET field of the file.",
    ).tag(config=True)

    pointing_dec = Float(
        default_value=None,
        allow_none=True,
        help="Pointing declination in deg. Overrides the TARGET field of the file.",
    ).tag(config=True)

    pointing_update_interval = Float(
        default_value=1.0,
        help=(
            "The altitude and azimuth of the pointing are recomputed when the time of the"
            " event differs by more than this value (in s) from the last computation."
        ),
    ).tag(config=True)

    focal_length_choice = UseEnum(
        FocalLengthKind,
        default_value=FocalLengthKind.EQUIVALENT,
        help="Which focal length to use for the camera frame transformations.",
    ).tag(config=True)

    swapped_modules = List(
        trait=Dict(),
        default_value=[],
        help=(
            "Wrongly connected modules: the pixels of the two modules are swapped for the"
            " events of the telescope in the period (UTC), e.g. [{\"tel_id\": 22,"
            " \"start\": \"2023-09-01T00:00:00\", \"stop\": \"2024-07-18T00:00:00\","
            " \"modules\": [59, 88]}]"
        ),
    ).tag(config=True)

    def __init__(self, input_url=None, config=None, parent=None, **kwargs):
        # LST/CTA uses differenct filename naming convention, how to work with the SST1M file naming convention?
        # A list of files can also be given as input_url, they are read one after the other.
        # ctapipe.io.EventSource only knows about the first one.
        input_urls = None
        if isinstance(input_url, list | tuple):
            input_urls = list(input_url)
            input_url = input_urls[0]

        super().__init__(input_url=input_url, config=config, parent=parent, **kwargs)

        self._input_urls = [self.input_url] if input_urls is None else [
            EventSource.input_url.validate(self, url) for url in input_urls
        ]
        # input files of the current provenance activity (zfits files have no reference metadata)
        for path in self.filelist:
            Provenance().add_input_file(path, role="R0/Event", add_meta=False)

        # obs_id from the date and run number of the file name
        date_run = parse_file_name(self.filelist[0])
        self.run_number = int(date_run[1]) if date_run else 0
        self.run_id = int(''.join(date_run)) if date_run else 0
        self.tel_id = 0

        # LST reads camera_config from input files, is it needed such functionality for SST1M?
        self.camera_config = None
        self.run_start = Time(self.camera_config.date, format='unix') if self.camera_config is not None else None

        self._subarray = SUBARRAY_DESCRIPTION


        # Target and pointing from the TARGET field of the file, unless given by the user
        header = fits.getheader(self.filelist[0], 'Events')
        self._target, self._wobble, ra, dec = parse_target_field(header.get('TARGET'))
        self._pointing_manual = (self.pointing_ra is not None) and (self.pointing_dec is not None)
        if self._pointing_manual:
            ra, dec = self.pointing_ra, self.pointing_dec

        self._pointing = None
        self._tel_locations = {}
        self._altaz_cache = {}
        target_info = {}
        pointing_mode = PointingMode.UNKNOWN
        if (ra is not None) and (dec is not None):
            self._pointing = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
            target_info["subarray_pointing_lon"] = ra * u.deg
            target_info["subarray_pointing_lat"] = dec * u.deg
            target_info["subarray_pointing_frame"] = CoordinateFrameType.ICRS
            pointing_mode = PointingMode.TRACK

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

        self._swat_event_ids_available = self.check_swat_event_ids_available(self.filelist)

        self._pixel_swaps = parse_swapped_modules(self.swapped_modules)

    @property
    def filelist(self):
        """All the files read by the source"""
        return [str(url) for url in self._input_urls]

    @property
    def subarray(self):
        return self._subarray

    @property
    def target(self):
        """Target name from the TARGET field of the file"""
        return self._target

    @property
    def wobble(self):
        """Wobble from the TARGET field of the file (``W<n>``, ``UNDEF`` or None)"""
        return self._wobble

    @property
    def pointing(self):
        """Pointing direction (ICRS) of the run, None if unknown"""
        return self._pointing

    @property
    def pointing_manual(self):
        """True if the pointing is given by the user and not read from the file"""
        return self._pointing_manual

    def _tel_location(self, tel_id):
        if tel_id not in self._tel_locations:
            locations = self.subarray.tel_coords.to_earth_location()
            self._tel_locations[tel_id] = locations[self.subarray.tel_index_array[tel_id]]
        return self._tel_locations[tel_id]

    def swapped_pixels(self, tel_id, local_camera_clock):
        """
        Pixel ids of the wrongly connected modules (pairs) of the telescope ``tel_id``
        at the time ``local_camera_clock`` (ns, TAI), see ``swapped_modules``
        """
        time = local_camera_clock / 1e9
        return [
            (pixels_1, pixels_2)
            for start, stop, pixels_1, pixels_2 in self._pixel_swaps.get(tel_id, [])
            if start < time < stop
        ]

    def _pixel_order(self, tel_id, pixel_ids, local_camera_clock):
        """
        Order of the pixels of the event in the camera: by pixel id, with the
        pixels of the wrongly connected modules swapped
        """
        order = np.argsort(pixel_ids)
        for pixels_1, pixels_2 in self.swapped_pixels(tel_id, local_camera_clock):
            order[pixels_1], order[pixels_2] = order[pixels_2], order[pixels_1]
        return order

    def _fill_trigger_and_pointing(self, array_event, tel_id, local_camera_clock):
        time = camera_clock_to_time(local_camera_clock)
        array_event.trigger.time = time
        array_event.trigger.tel[tel_id].time = time
        array_event.trigger.tels_with_trigger = [tel_id]

        if not self.pointing_information or self._pointing is None:
            return

        # the alt/az transformation is slow, it is only recomputed when the time changed enough
        cached = self._altaz_cache.get(tel_id)
        if cached is None or abs((time - cached[0]).to_value(u.s)) > self.pointing_update_interval:
            horizon_frame = AltAz(obstime=time, location=self._tel_location(tel_id))
            altaz = self._pointing.transform_to(horizon_frame)
            cached = (time, altaz.az.to(u.rad), altaz.alt.to(u.rad))
            self._altaz_cache[tel_id] = cached
        _, azimuth, altitude = cached

        pointing = array_event.pointing
        pointing.tel[tel_id].azimuth = azimuth
        pointing.tel[tel_id].altitude = altitude
        pointing.array_azimuth = azimuth
        pointing.array_altitude = altitude
        pointing.array_ra = self._pointing.ra.to(u.rad)
        pointing.array_dec = self._pointing.dec.to(u.rad)

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
        return (DataLevel.R0, )

    @property
    def swat_event_ids_available(self):
        return self._swat_event_ids_available

    @staticmethod
    def check_swat_event_ids_available(filelist, n_events=N_EVENTS_SWAT_ID_CHECK):
        """
        Determine if the files contain the array event ids (``arrayEvtNum``)
        written by SWAT.

        If SWAT did not write them, ``arrayEvtNum`` is always 0. Otherwise it can
        be 0 at most once, if SWAT was just restarted. The ids are thus considered
        available in a file if any of its first ``n_events`` events has a non zero
        ``arrayEvtNum``.

        Parameters
        ----------
        filelist: list of str or str
            Files of the run
        n_events: int
            Number of events read at the beginning of each file

        Returns
        -------
        bool:
            True if all the files contain the SWAT ids. If only some of them do,
            False is returned (with a warning) so that the event ids of the run
            are consistent.
        """
        if isinstance(filelist, str | os.PathLike):
            filelist = [filelist]

        available = [file_has_swat_event_ids(path, n_events) for path in filelist]

        if any(available) and not all(available):
            logger.warning(
                "SWAT event ids are available only in some of the files, they are not used: %s",
                {str(path): has_ids for path, has_ids in zip(filelist, available, strict=True)},
            )

        return len(available) > 0 and all(available)

    def _generator(self):
        """
        Read the files one after the other.
        NOTE: protozfits.MultiZFitsFiles merges interleaved files by event_id (LST),
        SST-1M files are written one after the other and have no event_id field
        """
        count = 0
        for input_path in self.filelist:
            for array_event in self.get_array_event(input_path):
                array_event.count = count
                array_event.index.obs_id = self.run_id

                yield array_event
                count += 1

    def get_array_event(self, input_path):
        """
        Read the events of a single file. Only the R0 data (``event.r0``),
        the trigger and the pointing are filled.
        """
        self.log.info("Reading %s", input_path)
        array_event = SST1MArrayEventContainer()
        with File(input_path) as f:
            array_event.r0.meta = dict(is_simulation=False)
            for event_counter, event in enumerate(f.Events):
                if self._swat_event_ids_available:
                    array_event.index.event_id = event.arrayEvtNum
                else:
                    array_event.index.event_id = event.eventNumber

                tel_id = event.telescopeID
                pixel_ids = event.hiGain.waveforms.pixelsIndices
                n_pixels = len(pixel_ids)
                local_camera_clock = (
                    np.int64(event.local_time_sec * 1E9) +
                    np.int64(event.local_time_nanosec)
                )
                sort_ids = self._pixel_order(tel_id, pixel_ids, local_camera_clock)
                samples = event.hiGain.waveforms.samples.reshape(n_pixels, -1)
                n_samples = samples.shape[1]

                try:
                    unsorted_baseline = event.hiGain.waveforms.baselines
                except AttributeError as err:
                    raise AttributeError("Could not read `hiGain.waveforms.baselines`"
                        f"for event:{event_counter} (eventNumber {event.eventNumber})\n"
                        f"of file:{input_path}\n") from err

                array_event.r0.tel.clear()
                r0 = array_event.r0.tel[tel_id]
                r0.waveform = samples[sort_ids].reshape(1, n_pixels, n_samples)
                r0.num_samples = n_samples
                r0.pedestal = unsorted_baseline[sort_ids] / 16
                r0.camera_event_number = event.eventNumber
                r0.pixel_flags = event.pixels_flags[sort_ids]
                r0.local_camera_clock = local_camera_clock
                if event.trig is not None:
                    r0.gps_time = (
                        np.int64(event.trig.timeSec * 1E9) +
                        np.int64(event.trig.timeNanoSec)
                    )
                else:
                    r0.gps_time = np.int64(0)
                r0.camera_event_type = event.event_type
                r0.array_event_type = event.eventType
                r0.trigger_input_traces = self._read_trigger_traces(
                    event.trigger_input_traces, self._prepare_trigger_input,
                    "trigger_input_traces", n_samples,
                )
                r0.trigger_output_patch7 = self._read_trigger_traces(
                    event.trigger_output_patch7, self._prepare_trigger_output,
                    "trigger_output_patch7", n_samples,
                )
                r0.trigger_output_patch19 = self._read_trigger_traces(
                    event.trigger_output_patch19, self._prepare_trigger_output,
                    "trigger_output_patch19", n_samples,
                )

                self._fill_trigger_and_pointing(array_event, tel_id, r0.local_camera_clock)
                # internal triggers are the pedestal events
                array_event.trigger.event_type = (
                    EventType.SKY_PEDESTAL if r0.camera_event_type == CameraEventType.INTERNAL
                    else EventType.SUBARRAY
                )
                yield array_event

    @staticmethod
    def _read_trigger_traces(traces, prepare, name, n_samples):
        if len(traces) > 0:
            return prepare(traces)
        warnings.warn(f'{name} does not exist: --> nan', stacklevel=3)
        return np.full((432, n_samples), np.nan)

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
        """
        SST-1M zfits files have an ``Events`` table of DigiCam protobuf messages
        """
        try:
            with fits.open(file_path) as hdul:
                if "Events" not in hdul:
                    return False
                return hdul["Events"].header.get("PBFHEAD") == "DataModel.CameraEvent"
        except (OSError, TypeError, ValueError):
            return False
