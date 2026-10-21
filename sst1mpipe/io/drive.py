"""
Logs of the drive system of the SST-1M telescopes (files ``DriveSystem<n>_<date>_<nnn>.fits``):
position of the telescope measured by the drive (altitude, azimuth) and state of the drive.
"""
import logging
import re

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, SkyCoord, frame_transform_graph
from astropy.io import fits
from astropy.table import QTable
from astropy.time import Time
from ctapipe.instrument import SubarrayDescription

from sst1mpipe.resources import SUBARRAY_FILE

__all__ = ["read_drive_logs", "drive_radec"]

logger = logging.getLogger(__name__)

# the drive system <n> (SUBSYS = DriveSystem<n>) is the one of the telescope 20 + n
DRIVE_TEL_ID_OFFSET = 20
# state of the drive in the logs
STATE_COLUMNS = ("is_tracking", "is_on_source", "is_moving", "is_in_park_position")


def read_drive_logs(file, subarray=None):
    """
    Position of the telescope measured by its drive system (``current_position_el`` and
    ``current_position_az``, deg) at the time of the drive (``current_time``, ms since 1970,
    in the time system ``TIMESSYS`` of the logs), and state of the drive (`STATE_COLUMNS`).
    The drive position includes the pointing model of the drive: the direction pointed to in the
    sky is given by `sst1mpipe.instrument.pointing_model.PointingModel.revert`.

    Parameters
    ----------
    file: str or Path
        Drive logs (``logs`` table)
    subarray: ctapipe.instrument.SubarrayDescription
        Location of the telescopes, by default the SST-1M subarray (``SUBARRAY_FILE``)

    Returns
    -------
    astropy.table.QTable
        ``time``, ``pointing`` (AltAz of the telescope location) and the state columns, with
        ``tel_id``, ``frame`` and ``equinox`` (sky frame of the drive) in the meta
    """
    if subarray is None:
        subarray = SubarrayDescription.from_hdf(SUBARRAY_FILE)

    with fits.open(file) as hdul:
        header = hdul["logs"].header
        data = hdul["logs"].data

        tel_id = int(re.search(r"\d+$", header["SUBSYS"]).group()) + DRIVE_TEL_ID_OFFSET
        frame = header["RADECSYS"].lower()
        equinox = f"J{header['EQUINOX']:.0f}"
        time_scale = header["TIMESSYS"].lower()
        t_start = Time(header["TSTART"], format="mjd", scale=header["TIMESYS"].lower())
        t_stop = Time(header["TSTOP"], format="mjd", scale=header["TIMESYS"].lower())
        logger.info("Drive logs of telescope %d from %s to %s", tel_id, t_start.isot, t_stop.isot)

        time = Time(data["current_time"] / 1000, format="unix", scale=time_scale)
        pointing = SkyCoord(
            alt=data["current_position_el"].astype(np.float64) * u.deg,
            az=data["current_position_az"].astype(np.float64) * u.deg,
            frame=AltAz(obstime=time, location=subarray.tel_earth_locations[tel_id]),
        )
        table = QTable(dict(time=time, pointing=pointing))
        for column in STATE_COLUMNS:
            table[column] = data[column].astype(bool)

    table.meta.update(tel_id=tel_id, frame=frame, equinox=equinox)
    return table


def drive_radec(drive):
    """
    Right ascension and declination of the position of the telescope of the drive logs
    (`read_drive_logs`), in the sky frame of the drive (``frame`` and ``equinox`` of the meta)
    """
    frame = frame_transform_graph.lookup_name(drive.meta["frame"])
    if "equinox" in frame.frame_attributes:
        frame = frame(equinox=drive.meta["equinox"])
    return drive["pointing"].transform_to(frame)
