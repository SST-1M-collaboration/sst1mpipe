"""
DL3 of the data (GADF event lists of gammapy) from the DL2 of ctapipe (``sst1mpipe-apply-models``),
with the cuts of ``ctapipe-optimize-event-selection`` and the point-like IRFs of
``ctapipe-compute-irf``, see ``sst1mpipe-dl2-dl3``.

One DL3 file per observation (``obs_id``): the events passing the cuts (``EVENTS``: reconstructed
energy, direction, gammaness), the good time interval (``GTI``), the pointing (``POINTING``), and the
IRFs (``EFFECTIVE AREA``, ``ENERGY DISPERSION``, ``RAD_MAX``, ``BACKGROUND``). The index files of
gammapy (``hdu-index.fits.gz``, ``obs-index.fits.gz``) are made by `write_index`.
"""
import operator
import re
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, SkyCoord
from astropy.io import fits
from astropy.table import QTable, Table
from astropy.time import Time
from ctapipe.core import QualityQuery
from pyirf.cuts import evaluate_binned_cut

__all__ = [
    "IRF_HDUS",
    "TIME_REFERENCE",
    "event_selection",
    "events_hdu",
    "gti_hdu",
    "irf_hdus",
    "pointing_hdu",
    "select_irf_directory",
    "write_index",
]

#: HDUs of the point-like IRFs of ctapipe-compute-irf copied to the DL3 files
IRF_HDUS = ["EFFECTIVE AREA", "ENERGY DISPERSION", "RAD_MAX", "BACKGROUND"]
#: reference of the times of the DL3 files (MJDREF, TT)
TIME_REFERENCE = Time(51910, format="mjd", scale="tt")


def select_irf_directory(directory, zenith):
    """
    IRF directory (``cuts.fits``, ``irf_point_like.fits.gz``) of the zenith angle closest to
    ``zenith``: ``directory`` itself if it has the IRFs, else its subdirectory ``ze<zenith in deg>``
    (e.g. ``irf/optimized/stereo``: ``ze20``, ``ze40``, ...) closest to ``zenith``
    """
    directory = Path(directory)
    if (directory / "irf_point_like.fits.gz").is_file():
        return directory
    zeniths = {}
    for path in directory.glob("ze*"):
        match = re.fullmatch(r"ze(\d+(?:\.\d+)?)", path.name)
        if match and (path / "irf_point_like.fits.gz").is_file():
            zeniths[float(match[1])] = path
    if not zeniths:
        raise FileNotFoundError(f"No IRFs in {directory} (irf_point_like.fits.gz, or ze<zenith>/irf_point_like.fits.gz)")
    zenith = u.Quantity(zenith, u.deg).to_value(u.deg)
    return zeniths[min(zeniths, key=lambda z: abs(z - zenith))]


def event_selection(events, cuts, energy_column, gammaness_column):
    """
    Mask of the events passing the cuts of ``ctapipe-optimize-event-selection``
    (`ctapipe.irf.OptimizationResult`): quality criteria, valid energy range, gammaness cut of their
    reconstructed energy bin (the direction cut, ``RAD_MAX``, is applied by the analysis)
    """
    criteria = [criterion for criterion in cuts.quality_query.quality_criteria if criterion != (" ", " ")]
    mask = QualityQuery(quality_criteria=criteria).get_table_mask(events) if criteria else np.ones(len(events), bool)
    energy = u.Quantity(events[energy_column], u.TeV)
    mask &= np.isfinite(energy) & np.isfinite(np.asarray(events[gammaness_column], dtype=float))
    mask &= (energy >= u.Quantity(cuts.valid_energy.min, u.TeV).min()) & (energy <= u.Quantity(cuts.valid_energy.max, u.TeV).max())
    gh = evaluate_binned_cut(np.asarray(events[gammaness_column], dtype=float), energy, cuts.gh_cuts, operator.ge)
    return mask & gh


def _time(time):
    """Seconds since `TIME_REFERENCE`"""
    return (Time(time).tt - TIME_REFERENCE).to_value(u.s)


def _time_header(header):
    header["MJDREFI"] = (int(TIME_REFERENCE.mjd), "[d] Integer part of the time reference MJD")
    header["MJDREFF"] = (TIME_REFERENCE.mjd % 1, "[d] Fractional part of the time reference MJD")
    header["TIMEUNIT"] = ("s", "Time unit")
    header["TIMESYS"] = ("TT", "Time system")
    header["TIMEREF"] = ("TOPOCENTER", "Time reference")


def events_hdu(events, obs_id, location, pointing, start, stop, meta):
    """
    ``EVENTS`` HDU of one observation

    Parameters
    ----------
    events: table with ``event_id``, ``time``, ``energy`` (TeV), ``alt``, ``az``, ``gammaness``
    obs_id: int
    location: astropy.coordinates.EarthLocation, of the array
    pointing: SkyCoord (ICRS) of the pointing
    start, stop: astropy.time.Time, of the good time interval
    meta: dict of header keywords (e.g. ``OBJECT``, ``TELESCOP``, ``ALT_PNT``, ``AZ_PNT``)
    """
    time = Time(events["time"])
    altaz = SkyCoord(alt=u.Quantity(events["alt"], u.deg), az=u.Quantity(events["az"], u.deg),
                     frame=AltAz(obstime=time, location=location))
    icrs = altaz.icrs
    table = Table()
    table["EVENT_ID"] = np.asarray(events["event_id"], dtype=np.int64)
    table["TIME"] = u.Quantity(_time(time), u.s)
    table["RA"] = icrs.ra.to(u.deg)
    table["DEC"] = icrs.dec.to(u.deg)
    table["ENERGY"] = u.Quantity(events["energy"], u.TeV)
    table["ALT"] = altaz.alt.to(u.deg)
    table["AZ"] = altaz.az.to(u.deg)
    table["GAMMANESS"] = np.asarray(events["gammaness"], dtype=np.float32)
    hdu = fits.BinTableHDU(table, name="EVENTS")
    header = hdu.header
    for key, value in (("HDUCLASS", "GADF"), ("HDUDOC", "https://gamma-astro-data-formats.readthedocs.io"),
                       ("HDUVERS", "0.3"), ("HDUCLAS1", "EVENTS"), ("OBS_ID", int(obs_id)),
                       ("OBS_MODE", "POINTING")):
        header[key] = value
    _time_header(header)
    livetime = (stop - start).to_value(u.s)
    header["TSTART"] = (_time(start), "[s] Start of the observation")
    header["TSTOP"] = (_time(stop), "[s] End of the observation")
    for key, time in (("OBS", start), ("END", stop)):
        date, clock = time.utc.isot.split("T")
        header[f"DATE-{key}"] = (date, "UTC date")
        header[f"TIME-{key}"] = (clock, "UTC time")
    header["ONTIME"] = (livetime, "[s] Total good time")
    # TODO: dead time of the cameras
    header["DEADC"] = (1.0, "Dead time correction factor (no dead time)")
    header["LIVETIME"] = (livetime, "[s] ONTIME * DEADC")
    header["RA_PNT"] = (pointing.ra.to_value(u.deg), "[deg] Pointing right ascension")
    header["DEC_PNT"] = (pointing.dec.to_value(u.deg), "[deg] Pointing declination")
    header["RADECSYS"] = ("ICRS", "Equatorial coordinate system")
    header["EQUINOX"] = (2000.0, "Equinox of the coordinates")
    geodetic = location.geodetic
    header["GEOLON"] = (geodetic.lon.to_value(u.deg), "[deg] Longitude of the array")
    header["GEOLAT"] = (geodetic.lat.to_value(u.deg), "[deg] Latitude of the array")
    header["ALTITUDE"] = (geodetic.height.to_value(u.m), "[m] Altitude of the array")
    header["EUNIT"] = "TeV"
    for key, value in meta.items():
        header[key] = value
    return hdu


def gti_hdu(start, stop):
    """``GTI`` HDU: one interval"""
    table = Table({"START": u.Quantity([_time(start)], u.s), "STOP": u.Quantity([_time(stop)], u.s)})
    hdu = fits.BinTableHDU(table, name="GTI")
    hdu.header["HDUCLASS"] = "GADF"
    hdu.header["HDUCLAS1"] = "GTI"
    _time_header(hdu.header)
    return hdu


def pointing_hdu(time, alt, az, location):
    """``POINTING`` HDU: the pointing of the telescopes (alt/az, RA/Dec) at the times ``time``"""
    time = Time(time)
    altaz = SkyCoord(alt=u.Quantity(alt, u.deg), az=u.Quantity(az, u.deg), frame=AltAz(obstime=time, location=location))
    table = QTable({"TIME": u.Quantity(_time(time), u.s), "RA_PNT": altaz.icrs.ra.to(u.deg),
                    "DEC_PNT": altaz.icrs.dec.to(u.deg), "ALT_PNT": altaz.alt.to(u.deg), "AZ_PNT": altaz.az.to(u.deg)})
    hdu = fits.BinTableHDU(table, name="POINTING")
    hdu.header["HDUCLASS"] = "GADF"
    hdu.header["HDUCLAS1"] = "POINTING"
    _time_header(hdu.header)
    return hdu


def irf_hdus(irf_file):
    """The HDUs `IRF_HDUS` of a file of ``ctapipe-compute-irf`` (those found)"""
    with fits.open(irf_file) as hdul:
        return [hdul[name].copy() for name in IRF_HDUS if name in hdul]


#: HDU_TYPE of the HDU index of gammapy per HDUCLAS2 of the IRF HDUs (GADF)
_IRF_TYPES = {"EFF_AREA": "aeff", "EDISP": "edisp", "RAD_MAX": "rad_max", "BKG": "bkg", "PSF": "psf"}


def _hdu_rows(path, directory):
    """Rows of the HDU index of gammapy of the HDUs of a DL3 file"""
    file_dir = str(path.parent.absolute().relative_to(directory.absolute(), walk_up=True))
    rows = []
    with fits.open(path) as hdul:
        obs_id = hdul["EVENTS"].header["OBS_ID"]
        for hdu in hdul[1:]:
            header = hdu.header
            if hdu.name in ("EVENTS", "GTI", "POINTING"):
                hdu_type = hdu_class = hdu.name.lower()
            elif header.get("HDUCLAS2") in _IRF_TYPES:
                hdu_type, hdu_class = _IRF_TYPES[header["HDUCLAS2"]], header.get("HDUCLAS4", "").lower()
            else:
                continue
            rows.append({"OBS_ID": obs_id, "HDU_TYPE": hdu_type, "HDU_CLASS": hdu_class, "FILE_DIR": file_dir,
                         "FILE_NAME": path.name, "HDU_NAME": hdu.name})
    return rows


def write_index(dl3_files, directory):
    """
    Index files of gammapy (``hdu-index.fits.gz``, ``obs-index.fits.gz``) of DL3 files with their
    IRFs (point-like IRFs: ``RAD_MAX``, 2D background), in ``directory`` (paths relative to it)

    Returns
    -------
    gammapy.data.DataStore
    """
    from gammapy.data import DataStore
    from gammapy.data.hdu_index_table import HDUIndexTable

    dl3_files = [Path(path) for path in dl3_files]
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    hdu_table = HDUIndexTable(rows=[row for path in dl3_files for row in _hdu_rows(path, directory)],
                              names=["OBS_ID", "HDU_TYPE", "HDU_CLASS", "FILE_DIR", "FILE_NAME", "HDU_NAME"])
    hdu_table.meta.update({"HDUCLASS": "GADF", "HDUDOC": "https://gamma-astro-data-formats.readthedocs.io",
                           "HDUVERS": "0.3", "HDUCLAS1": "INDEX", "HDUCLAS2": "HDU"})
    hdu_table.write(directory / "hdu-index.fits.gz", overwrite=True)
    # observation index of gammapy, from the EVENTS headers
    obs_table = DataStore.from_events_files([path.absolute() for path in dl3_files],
                                            irfs_paths=[path.absolute() for path in dl3_files]).obs_table
    obs_table.remove_columns([name for name in ("EVENTS_FILENAME", "IRF_FILENAME") if name in obs_table.colnames])
    pointing = {}
    for path in dl3_files:
        header = fits.getheader(path, "EVENTS")
        pointing[header["OBS_ID"]] = (header["ALT_PNT"], header["AZ_PNT"])
    obs_table["ALT_PNT"] = [pointing[obs_id][0] for obs_id in obs_table["OBS_ID"]] * u.deg
    obs_table["AZ_PNT"] = [pointing[obs_id][1] for obs_id in obs_table["OBS_ID"]] * u.deg
    obs_table["ZEN_PNT"] = 90 * u.deg - obs_table["ALT_PNT"].quantity
    obs_table.write(directory / "obs-index.fits.gz", overwrite=True)
    return DataStore.from_dir(directory)
