"""
DL3 of the data in the Gamma Astro Data Formats (GADF, https://gamma-astro-data-formats.readthedocs.io,
read by gammapy) from the DL2 of ctapipe (``sst1mpipe-apply-models``),
with the cuts of ``ctapipe-optimize-event-selection`` and the point-like IRFs of
``ctapipe-compute-irf``, see ``sst1mpipe-dl2-dl3``.

One DL3 file per observation (``obs_id``): the events passing the cuts (``EVENTS``: reconstructed
energy, direction, gammaness), the good time interval (``GTI``), the pointing (``POINTING``), and the
IRFs (``EFFECTIVE AREA``, ``ENERGY DISPERSION``, ``RAD_MAX``, ``BACKGROUND``). The index files
(``hdu-index.fits.gz``: ``HDU_INDEX``, ``obs-index.fits.gz``: ``OBS_INDEX``) are made by `write_index`.
`check_dl3_file` and `check_index` list what a DL3 file or index files miss of the GADF specification.
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
    "GADF_HDUDOC",
    "GADF_HDUVERS",
    "IRF_HDUS",
    "check_dl3_file",
    "check_index",
    "TIME_REFERENCE",
    "event_selection",
    "events_hdu",
    "gti_hdu",
    "irf_hdus",
    "obs_index_table",
    "pointing_hdu",
    "select_irf_directory",
    "write_index",
]

#: HDUs of the point-like IRFs of ctapipe-compute-irf copied to the DL3 files
IRF_HDUS = ["EFFECTIVE AREA", "ENERGY DISPERSION", "RAD_MAX", "BACKGROUND"]
#: reference of the times of the DL3 files (MJDREF, TT)
TIME_REFERENCE = Time(51910, format="mjd", scale="tt")
#: HDUDOC and HDUVERS of the GADF specification
GADF_HDUDOC = "https://github.com/open-gamma-ray-astro/gamma-astro-data-formats"
GADF_HDUVERS = "0.3"


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
    """Reference time of the TIME / START / STOP / TSTART / TSTOP values (GADF general/time)"""
    header["MJDREFI"] = (int(TIME_REFERENCE.mjd), "[d] Integer part of the time reference MJD")
    header["MJDREFF"] = (TIME_REFERENCE.mjd % 1, "[d] Fractional part of the time reference MJD")
    header["TIMEUNIT"] = ("s", "Time unit")
    header["TIMESYS"] = ("TT", "Time system")
    header["TIMEREF"] = ("TOPOCENTER", "Time reference")


def _gadf_header(header, hduclas1, hduclas2=None):
    """HDUCLASS, HDUDOC, HDUVERS, HDUCLAS1 (and HDUCLAS2) of the GADF specification"""
    header["HDUCLASS"] = ("GADF", "Gamma Astro Data Formats")
    header["HDUDOC"] = GADF_HDUDOC  # no comment: the URL fills the card
    header["HDUVERS"] = (GADF_HDUVERS, "Version of the format")
    header["HDUCLAS1"] = (hduclas1, "Type of HDU")
    if hduclas2 is not None:
        header["HDUCLAS2"] = (hduclas2, "Type of HDU")


def _location_header(header, location):
    """Earth location of the observatory (GEOLON, GEOLAT, ALTITUDE, as gammapy)"""
    geodetic = location.geodetic
    header["GEOLON"] = (geodetic.lon.to_value(u.deg), "[deg] Longitude of the array")
    header["GEOLAT"] = (geodetic.lat.to_value(u.deg), "[deg] Latitude of the array")
    header["ALTITUDE"] = (geodetic.height.to_value(u.m), "[m] Altitude of the array")


def _creator():
    import sst1mpipe

    return f"sst1mpipe {sst1mpipe.__version__}"


def events_hdu(events, obs_id, location, pointing, start, stop, meta):
    """
    ``EVENTS`` HDU of one observation (GADF events/events): the mandatory columns ``EVENT_ID``,
    ``TIME``, ``RA``, ``DEC``, ``ENERGY`` and the optional ``ALT``, ``AZ``, ``GAMMANESS`` (and ``MULTIP``),
    the mandatory header keywords (``ORIGIN``, ``TELESCOP``, ``INSTRUME``, ``CREATOR`` from ``meta``, with
    defaults)

    Parameters
    ----------
    events: table with ``event_id``, ``time``, ``energy`` (TeV), ``alt``, ``az``, ``gammaness``, and
        optionally ``multiplicity`` (number of telescopes of the reconstruction)
    obs_id: int
    location: astropy.coordinates.EarthLocation, of the array
    pointing: SkyCoord (ICRS) of the pointing
    start, stop: astropy.time.Time, of the good time interval
    meta: dict of header keywords (e.g. ``OBJECT``, ``RA_OBJ``, ``DEC_OBJ``, ``ORIGIN``, ``TELESCOP``,
        ``INSTRUME``, ``ALT_PNT``, ``AZ_PNT``)
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
    if "multiplicity" in events:
        table["MULTIP"] = np.asarray(events["multiplicity"], dtype=np.int16)
    hdu = fits.BinTableHDU(table, name="EVENTS")
    header = hdu.header
    _gadf_header(header, "EVENTS")
    header["OBS_ID"] = (int(obs_id), "Observation identifier")
    header["OBS_MODE"] = ("POINTING", "Observation mode (wobble: pointing)")
    _time_header(header)
    livetime = (stop - start).to_value(u.s)
    header["TSTART"] = (_time(start), "[s] Start of the observation")
    header["TSTOP"] = (_time(stop), "[s] End of the observation")
    header["TELAPSE"] = (_time(stop) - _time(start), "[s] TSTOP - TSTART")
    # human readable times, in the time system TIMESYS
    for key, time in (("OBS", start), ("END", stop)):
        date, clock = time.tt.isot.split("T")
        header[f"DATE-{key}"] = (date, "Date (TT)")
        header[f"TIME-{key}"] = (clock, "Time (TT)")
    header["DATE-BEG"] = (start.tt.isot, "Start of the observation (TT)")
    header["ONTIME"] = (livetime, "[s] Total good time")
    # TODO: dead time of the cameras, GTI of the observation (issue #252)
    header["DEADC"] = (1.0, "Dead time correction factor (no dead time)")
    header["LIVETIME"] = (livetime, "[s] ONTIME * DEADC")
    header["RA_PNT"] = (pointing.ra.to_value(u.deg), "[deg] Pointing right ascension")
    header["DEC_PNT"] = (pointing.dec.to_value(u.deg), "[deg] Pointing declination")
    header["RADECSYS"] = ("ICRS", "Equatorial coordinate system")
    header["EQUINOX"] = (2000.0, "Equinox of the coordinates")
    _location_header(header, location)
    header["EUNIT"] = "TeV"
    header["ORIGIN"] = ("SST-1M Collaboration", "Organisation that created the file")
    header["TELESCOP"] = ("SST-1M", "Telescope")
    header["INSTRUME"] = ("SST-1M", "Instrument")
    header["CREATOR"] = (_creator(), "Software that created the file")
    header["CREATED"] = (Time.now().utc.isot.split(".")[0], "Creation date (UTC)")
    for key, value in meta.items():
        header[key] = value
    return hdu


def gti_hdu(start, stop):
    """``GTI`` HDU: one interval"""
    table = Table({"START": u.Quantity([_time(start)], u.s), "STOP": u.Quantity([_time(stop)], u.s)})
    hdu = fits.BinTableHDU(table, name="GTI")
    _gadf_header(hdu.header, "GTI")
    _time_header(hdu.header)
    return hdu


def pointing_hdu(time, alt, az, location):
    """``POINTING`` HDU: the pointing of the telescopes (alt/az, RA/Dec) at the times ``time``"""
    time = Time(time)
    altaz = SkyCoord(alt=u.Quantity(alt, u.deg), az=u.Quantity(az, u.deg), frame=AltAz(obstime=time, location=location))
    table = QTable({"TIME": u.Quantity(_time(time), u.s), "RA_PNT": altaz.icrs.ra.to(u.deg),
                    "DEC_PNT": altaz.icrs.dec.to(u.deg), "ALT_PNT": altaz.alt.to(u.deg), "AZ_PNT": altaz.az.to(u.deg)})
    hdu = fits.BinTableHDU(table, name="POINTING")
    _gadf_header(hdu.header, "POINTING")
    _time_header(hdu.header)
    _location_header(hdu.header, location)
    return hdu


def irf_hdus(irf_file):
    """The HDUs `IRF_HDUS` of a file of ``ctapipe-compute-irf`` (those found)"""
    with fits.open(irf_file) as hdul:
        return [hdul[name].copy() for name in IRF_HDUS if name in hdul]


#: HDU_TYPE of the HDU index per HDUCLAS2 of the IRF HDUs (GADF data_storage/hdu_index)
_IRF_TYPES = {"EFF_AREA": "aeff", "EDISP": "edisp", "RAD_MAX": "rad_max", "BKG": "bkg", "PSF": "psf", "RPSF": "psf"}
#: HDU_TYPE and HDU_CLASS of the GADF specification
GADF_HDU_TYPES = {"events", "gti", "aeff", "psf", "edisp", "bkg", "rad_max"}
GADF_HDU_CLASSES = {"events", "gti", "aeff_2d", "aeff_3d", "edisp_2d", "psf_table", "psf_3gauss", "psf_king",
                    "bkg_2d", "bkg_3d", "rad_max_2d"}


def _hdu_rows(path, directory):
    """
    Rows of the HDU index of the HDUs of a DL3 file: EVENTS, GTI and the IRFs (the HDU types of the
    GADF specification; the POINTING HDU is in the files, not in the index)
    """
    file_dir = str(path.parent.absolute().relative_to(directory.absolute(), walk_up=True))
    rows = []
    with fits.open(path) as hdul:
        obs_id = hdul["EVENTS"].header["OBS_ID"]
        for hdu in hdul[1:]:
            header = hdu.header
            if hdu.name in ("EVENTS", "GTI"):
                hdu_type = hdu_class = hdu.name.lower()
            elif header.get("HDUCLAS1") == "RESPONSE" and header.get("HDUCLAS2") in _IRF_TYPES:
                hdu_type, hdu_class = _IRF_TYPES[header["HDUCLAS2"]], header.get("HDUCLAS4", "").lower()
            else:
                continue
            rows.append({"OBS_ID": obs_id, "HDU_TYPE": hdu_type, "HDU_CLASS": hdu_class, "FILE_DIR": file_dir,
                         "FILE_NAME": path.name, "HDU_NAME": hdu.name, "SIZE": path.stat().st_size})
    return rows


#: columns of the observation index (GADF data_storage/obs_index) from the EVENTS headers: name, unit
_OBS_COLUMNS = [
    ("OBS_ID", None), ("OBS_MODE", None), ("RA_PNT", u.deg), ("DEC_PNT", u.deg), ("ALT_PNT", u.deg),
    ("AZ_PNT", u.deg), ("TSTART", u.s), ("TSTOP", u.s), ("ONTIME", u.s), ("LIVETIME", u.s), ("DEADC", None),
    ("DATE-OBS", None), ("TIME-OBS", None), ("DATE-END", None), ("TIME-END", None), ("N_TELS", None),
    ("TELLIST", None), ("OBJECT", None), ("RA_OBJ", u.deg), ("DEC_OBJ", u.deg),
]


def obs_index_table(dl3_files):
    """
    Observation index (GADF data_storage/obs_index, HDU ``OBS_INDEX``) of DL3 files, from their
    ``EVENTS`` headers: the required columns ``OBS_ID``, ``OBS_MODE``, ``RA_PNT``, ``DEC_PNT``,
    ``TSTART``, ``TSTOP`` and optional ones (``ALT_PNT``, ``AZ_PNT``, ``ZEN_PNT``, ``ONTIME``,
    ``LIVETIME``, ``DEADC``, dates, telescopes, ``OBJECT``, ``TMID``, ``EVENT_COUNT``, ``QUALITY``);
    header: GADF keywords, reference time and location of the observatory (of the first file, the same
    for all)
    """
    headers = [fits.getheader(path, "EVENTS") for path in dl3_files]
    rows = []
    for header in headers:
        for key in ("MJDREFI", "MJDREFF", "TIMESYS"):
            if header.get(key) != headers[0].get(key):
                raise ValueError(f"OBS_ID {header['OBS_ID']}: {key} differs from the other DL3 files")
        rows.append(header)
    table = Table()
    for name, unit in _OBS_COLUMNS:
        if not all(name in header for header in rows):
            continue
        values = [header[name] for header in rows]
        table[name] = u.Quantity(values, unit) if unit is not None else values
    if "ALT_PNT" in table.colnames:
        table["ZEN_PNT"] = 90 * u.deg - table["ALT_PNT"].quantity
    table["TMID"] = (table["TSTART"].quantity + table["TSTOP"].quantity) / 2
    table["EVENT_COUNT"] = [header["NAXIS2"] for header in rows]
    table["QUALITY"] = np.full(len(rows), -1, dtype=np.int16)  # unknown
    table.meta["EXTNAME"] = "OBS_INDEX"
    header = fits.Header()
    _gadf_header(header, "INDEX", "OBS")
    if rows:
        for key in ("MJDREFI", "MJDREFF", "TIMEUNIT", "TIMESYS", "TIMEREF", "GEOLON", "GEOLAT", "ALTITUDE",
                    "TELESCOP", "ORIGIN"):
            if key in rows[0]:
                header[key] = rows[0][key]
    table.meta.update(header)
    return table


def write_index(dl3_files, directory):
    """
    Index files (GADF data_storage) of DL3 files with their point-like IRFs (``RAD_MAX``, 2D
    background), in ``directory``, paths relative to it: ``hdu-index.fits.gz`` (``HDU_INDEX``) and
    ``obs-index.fits.gz`` (``OBS_INDEX``, `obs_index_table`)

    Returns
    -------
    gammapy.data.DataStore
    """
    from gammapy.data import DataStore

    dl3_files = sorted(Path(path) for path in dl3_files)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    hdu_table = Table(rows=[row for path in dl3_files for row in _hdu_rows(path, directory)] or None,
                      names=["OBS_ID", "HDU_TYPE", "HDU_CLASS", "FILE_DIR", "FILE_NAME", "HDU_NAME", "SIZE"],
                      dtype=[np.int64, str, str, str, str, str, np.int64])
    hdu_table["SIZE"].unit = u.byte
    hdu_table.meta["EXTNAME"] = "HDU_INDEX"
    header = fits.Header()
    _gadf_header(header, "INDEX", "HDU")
    hdu_table.meta.update(header)
    hdu_table.write(directory / "hdu-index.fits.gz", overwrite=True)
    obs_index_table(dl3_files).write(directory / "obs-index.fits.gz", overwrite=True)
    return DataStore.from_dir(directory)


#: mandatory columns and header keywords of the GADF specification per HDU (HDUCLAS1)
_MANDATORY = {
    "EVENTS": (["EVENT_ID", "TIME", "RA", "DEC", "ENERGY"],
               ["HDUCLASS", "HDUDOC", "HDUVERS", "HDUCLAS1", "OBS_ID", "TSTART", "TSTOP", "ONTIME", "LIVETIME",
                "DEADC", "OBS_MODE", "RA_PNT", "DEC_PNT", "EQUINOX", "RADECSYS", "ORIGIN", "TELESCOP", "INSTRUME",
                "CREATOR", "MJDREFI", "MJDREFF", "TIMEUNIT", "TIMESYS", "TIMEREF", "GEOLON", "GEOLAT", "ALTITUDE"]),
    "GTI": (["START", "STOP"], ["HDUCLASS", "HDUCLAS1", "MJDREFI", "MJDREFF", "TIMEUNIT", "TIMESYS", "TIMEREF"]),
    "POINTING": (["TIME", "RA_PNT", "DEC_PNT"],
                 ["HDUCLASS", "HDUCLAS1", "MJDREFI", "MJDREFF", "TIMEUNIT", "TIMESYS", "TIMEREF", "GEOLON",
                  "GEOLAT", "ALTITUDE"]),
}
_IRF_COLUMNS = {"AEFF_2D": ["ENERG_LO", "ENERG_HI", "THETA_LO", "THETA_HI", "EFFAREA"],
                "EDISP_2D": ["ENERG_LO", "ENERG_HI", "MIGRA_LO", "MIGRA_HI", "THETA_LO", "THETA_HI", "MATRIX"],
                "BKG_2D": ["ENERG_LO", "ENERG_HI", "THETA_LO", "THETA_HI", "BKG"],
                "RAD_MAX_2D": ["ENERG_LO", "ENERG_HI", "THETA_LO", "THETA_HI", "RAD_MAX"]}


def _check_hdu(hdu, columns, keywords):
    problems = [f"{hdu.name}: no column {name}" for name in columns if name not in hdu.columns.names]
    problems += [f"{hdu.name}: no keyword {key}" for key in keywords if key not in hdu.header]
    return problems


def check_dl3_file(path):
    """
    What a DL3 file misses of the GADF specification (mandatory columns and header keywords of
    ``EVENTS``, ``GTI``, ``POINTING``; point-like IRFs: ``HDUCLAS1-4``, columns, ``RAD_MAX_2D`` or the
    keyword ``RAD_MAX``; float64 times; ``DEADC`` in [0, 1])

    Returns
    -------
    list of str: the problems (empty: the file follows the specification)
    """
    problems = []
    with fits.open(path) as hdul:
        names = [hdu.name for hdu in hdul]
        for name in ("EVENTS", "GTI"):
            if name not in names:
                problems.append(f"no {name} HDU")
        for name, (columns, keywords) in _MANDATORY.items():
            if name in names:
                problems += _check_hdu(hdul[name], columns, keywords)
        for name, column in (("EVENTS", "TIME"), ("GTI", "START"), ("GTI", "STOP"), ("POINTING", "TIME")):
            if name in names and column in hdul[name].columns.names and hdul[name].columns[column].format != "D":
                problems.append(f"{name}: {column} is not float64")
        if "EVENTS" in names:
            header = hdul["EVENTS"].header
            if not 0 <= header.get("DEADC", -1) <= 1:
                problems.append("EVENTS: DEADC not in [0, 1]")
            if header.get("OBS_MODE") not in ("POINTING", "RASTER", "SLEW", "SCAN", "DRIFT"):
                problems.append(f"EVENTS: OBS_MODE {header.get('OBS_MODE')}")
        point_like, rad_max = False, False
        for hdu in hdul[1:]:
            header = hdu.header
            if header.get("HDUCLAS1") != "RESPONSE":
                continue
            problems += _check_hdu(hdu, _IRF_COLUMNS.get(header.get("HDUCLAS4"), []),
                                   ["HDUCLASS", "HDUDOC", "HDUVERS", "HDUCLAS2", "HDUCLAS3", "HDUCLAS4"])
            if header.get("HDUCLAS3") not in ("POINT-LIKE", "FULL-ENCLOSURE"):
                problems.append(f"{hdu.name}: HDUCLAS3 {header.get('HDUCLAS3')}")
            if header.get("HDUCLAS4", "").lower() not in GADF_HDU_CLASSES:
                problems.append(f"{hdu.name}: HDUCLAS4 {header.get('HDUCLAS4')}")
            point_like |= header.get("HDUCLAS3") == "POINT-LIKE"
            rad_max |= header.get("HDUCLAS4") == "RAD_MAX_2D" or "RAD_MAX" in header
        if point_like and not rad_max:
            problems.append("point-like IRFs without RAD_MAX_2D HDU nor RAD_MAX keyword")
    return problems


def check_index(directory):
    """
    What the index files of a directory miss of the GADF specification (``HDU_INDEX``: columns,
    ``HDU_TYPE`` and ``HDU_CLASS`` values, files; ``OBS_INDEX``: required columns; header keywords)

    Returns
    -------
    list of str: the problems (empty: the index files follow the specification)
    """
    directory = Path(directory)
    problems = []
    for file, extname, hduclas2, columns in (
        ("hdu-index.fits.gz", "HDU_INDEX", "HDU", ["OBS_ID", "HDU_TYPE", "HDU_CLASS", "FILE_DIR", "FILE_NAME", "HDU_NAME"]),
        ("obs-index.fits.gz", "OBS_INDEX", "OBS", ["OBS_ID", "OBS_MODE", "RA_PNT", "DEC_PNT", "TSTART", "TSTOP"]),
    ):
        path = directory / file
        if not path.is_file():
            problems.append(f"no {file}")
            continue
        with fits.open(path) as hdul:
            if extname not in hdul:
                problems.append(f"{file}: no HDU {extname}")
                continue
            hdu = hdul[extname]
            problems += _check_hdu(hdu, columns, ["HDUCLASS", "HDUDOC", "HDUVERS", "HDUCLAS1", "HDUCLAS2"])
            if hdu.header.get("HDUCLAS1") != "INDEX" or hdu.header.get("HDUCLAS2") != hduclas2:
                problems.append(f"{file}: HDUCLAS1/2 not INDEX/{hduclas2}")
            if extname == "HDU_INDEX":
                table = Table(hdu.data)
                for row in table:
                    if row["HDU_TYPE"] not in GADF_HDU_TYPES:
                        problems.append(f"{file}: HDU_TYPE {row['HDU_TYPE']}")
                    if row["HDU_CLASS"] not in GADF_HDU_CLASSES:
                        problems.append(f"{file}: HDU_CLASS {row['HDU_CLASS']}")
                    if not (directory / row["FILE_DIR"] / row["FILE_NAME"]).is_file():
                        problems.append(f"{file}: no file {row['FILE_DIR']}/{row['FILE_NAME']}")
    return problems
