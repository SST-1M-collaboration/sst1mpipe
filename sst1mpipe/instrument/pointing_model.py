"""
Pointing models of the SST-1M telescopes: corrections between the position of the drive
(encoders) and the direction pointed to in the sky, as a function of the altitude and azimuth.
The drive applies its pointing model to the position of the target, so its logs
(`sst1mpipe.io.drive.read_drive_logs`) give the drive position, and `PointingModel.revert` gives
the direction pointed to in the sky.
"""
import logging
from abc import abstractmethod

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import Dict, Float, FloatTelescopeParameter, List

__all__ = ["PointingModel", "TPointPointingModel"]

logger = logging.getLogger(__name__)


class PointingModel(TelescopeComponent):
    """
    Base class of the pointing models: corrections (d_alt, d_az) of the drive position of a
    telescope, sky position = drive position + corrections(drive position). A child class
    defines `corrections`.

    `revert` (drive -> sky) is exact. `apply` (sky -> drive) inverts it once, to first order
    (corrections taken at the sky position), with a warning if `revert` of the result is farther
    than ``tolerance`` from the sky position.
    """

    tolerance = Float(
        default_value=30.0,
        help="Largest distance (arcsec) between the sky position and revert of the drive position found by apply",
    ).tag(config=True)

    @abstractmethod
    def corrections(self, tel_id, alt, az, time=None):
        """
        Corrections (d_alt, d_az) of the drive position alt, az (Quantity, angles) of the
        telescope tel_id at ``time`` (astropy Time, for coefficients changing with time)
        """

    def revert(self, tel_id, pointing):
        """
        Sky position (AltAz SkyCoord) pointed to by the telescope tel_id at the drive position
        ``pointing`` (AltAz SkyCoord). Altitudes beyond +-90 deg are NaN.
        """
        return self._add_corrections(tel_id, pointing, 1)

    def apply(self, tel_id, pointing):
        """
        Drive position (AltAz SkyCoord) of the telescope tel_id pointing to the sky position
        ``pointing`` (AltAz SkyCoord): first order inverse of `revert`. A warning is logged if
        `revert` of the drive position is farther than ``tolerance`` from the sky position.
        Altitudes beyond +-90 deg are NaN.
        """
        drive = self._add_corrections(tel_id, pointing, -1)
        distance = np.atleast_1d(self.revert(tel_id, drive).separation(pointing).to_value(u.arcsec))
        above = distance > self.tolerance
        if np.any(above):
            logger.warning(
                "Pointing model of telescope %d: %d positions inverted with an error above %.1f arcsec (largest %.1f arcsec)",
                tel_id, np.count_nonzero(above), self.tolerance, np.max(distance[above]),
            )
        return drive

    def _add_corrections(self, tel_id, pointing, sign):
        """``pointing`` + sign * corrections at ``pointing``, at the obstime of its frame"""
        d_alt, d_az = self.corrections(tel_id, pointing.alt, pointing.az, time=pointing.frame.obstime)
        alt = (pointing.alt + sign * d_alt).to_value(u.deg)
        alt = np.where(np.abs(alt) <= 90, alt, np.nan)
        return SkyCoord(alt=alt * u.deg, az=pointing.az + sign * d_az, frame=pointing.frame)


class TPointPointingModel(PointingModel):
    """
    Pointing model of alt-az mounts with the terms of TPOINT [1]_ (names and units), corrections of
    the drive (encoder) position alt, az:

    - d_az = -IA - AN sin(az) tan(alt) - AW cos(az) tan(alt) - CA / cos(alt) - NPAE tan(alt)
      - ACES sin(az) + ACEC cos(az)
    - d_alt = IE - AN cos(az) + AW sin(az) - TF cos(alt) - TX / tan(alt) + ECES sin(alt)
      + ECEC cos(alt)

    The terms IA, IE, AN, AW, CA, NPAE, TF and TX have the signs of pytelpoint [2]_
    (``pytelpoint.transform.azel_model``). The centring terms ACES, ACEC, ECES and ECEC have the
    signs of the pointing model of the SST-1M drives (``AzElCorrCoefTel<n>.conf``, without TF and
    TX), in which d_alt has -IE: the drive positions of the night of 2026-01-20 need +IE.

    The coefficients are in arcsec, 0 by default (no correction), configured with the TPOINT names
    in lower case, e.g. ``{"TPointPointingModel": {"ia": [["id", 21, -8441.94], ["id", 22, -2074.65]]}}``.
    The coefficients loaded in a drive change with time: ``periods`` gives the coefficients of a
    telescope over periods of time, chosen with the time of the positions (obstime of their AltAz
    frame). A position of a telescope with periods outside all of them uses the coefficients of
    the parameters, with a warning. ``SST1M_POINTING_MODEL_CONFIG_FILE`` has the coefficients of
    the SST-1M drives.

    References
    ----------
    .. [1] P. T. Wallace, "TPOINT - Telescope Pointing Analysis System", Starlink User Note 100
       (1994)
    .. [2] T. E. Pickering, "pytelpoint: an open-source package for modeling and assessing
       telescope pointing performance", Proc. SPIE 12189, 121892L (2022),
       https://doi.org/10.1117/12.2630648
    """

    ia = FloatTelescopeParameter(default_value=0.0, help="IA: azimuth index (arcsec)").tag(config=True)
    ie = FloatTelescopeParameter(default_value=0.0, help="IE: elevation index (arcsec)").tag(config=True)
    an = FloatTelescopeParameter(
        default_value=0.0, help="AN: north-south misalignment of the azimuth axis (arcsec)"
    ).tag(config=True)
    aw = FloatTelescopeParameter(
        default_value=0.0, help="AW: east-west misalignment of the azimuth axis (arcsec)"
    ).tag(config=True)
    ca = FloatTelescopeParameter(
        default_value=0.0,
        help="CA: left-right collimation, non perpendicularity of the pointing direction and the elevation axis (arcsec)",
    ).tag(config=True)
    npae = FloatTelescopeParameter(
        default_value=0.0, help="NPAE: non perpendicularity of the azimuth and elevation axes (arcsec)"
    ).tag(config=True)
    tf = FloatTelescopeParameter(default_value=0.0, help="TF: tube flexure in cos(alt) (arcsec)").tag(config=True)
    tx = FloatTelescopeParameter(default_value=0.0, help="TX: tube flexure in 1 / tan(alt) (arcsec)").tag(config=True)
    aces = FloatTelescopeParameter(default_value=0.0, help="ACES: azimuth centring error, sin component (arcsec)").tag(config=True)
    acec = FloatTelescopeParameter(default_value=0.0, help="ACEC: azimuth centring error, cos component (arcsec)").tag(config=True)
    eces = FloatTelescopeParameter(default_value=0.0, help="ECES: elevation centring error, sin component (arcsec)").tag(config=True)
    ecec = FloatTelescopeParameter(default_value=0.0, help="ECEC: elevation centring error, cos component (arcsec)").tag(config=True)

    periods = List(
        trait=Dict(),
        default_value=[],
        help=(
            "Coefficients (arcsec) of a telescope over a period of time, instead of those of the parameters:"
            ' {"tel_id": 22, "start": "2024-04-30T12:00:00", "stop": "2026-03-14T12:00:00", "ia": -2103.0, ...}.'
            " start (included) and stop (excluded) are UTC, null for an open period. The missing terms are 0."
        ),
    ).tag(config=True)

    TERMS = ("ia", "ie", "an", "aw", "ca", "npae", "tf", "tx", "aces", "acec", "eces", "ecec")

    def __init__(self, subarray, config=None, parent=None, **kwargs):
        super().__init__(subarray, config=config, parent=parent, **kwargs)
        self._periods = {}
        for period in self.periods:
            unknown = set(period) - {"tel_id", "start", "stop", *self.TERMS}
            if "tel_id" not in period or unknown:
                raise ValueError(f"Pointing model period without tel_id or with unknown keys {sorted(unknown)}: {period}")
            start, stop = (Time(period[key], scale="utc") if period.get(key) is not None else None for key in ("start", "stop"))
            if start is not None and stop is not None and start >= stop:
                raise ValueError(f"Pointing model period starting after its stop: {period}")
            coefficients = {name: float(period.get(name, 0.0)) for name in self.TERMS}
            self._periods.setdefault(int(period["tel_id"]), []).append((start, stop, coefficients))
        for tel_id, periods in self._periods.items():
            periods.sort(key=lambda period: -np.inf if period[0] is None else period[0].mjd)
            for (_, stop, _), (start, _, _) in zip(periods[:-1], periods[1:], strict=True):
                if stop is None or start is None or stop > start:
                    raise ValueError(f"Overlapping pointing model periods of telescope {tel_id}")

    def coefficients(self, tel_id, time=None):
        """
        Coefficients (arcsec) of the telescope tel_id at ``time`` (astropy Time, scalar or array):
        those of the period including it, otherwise those of the parameters (with a warning if the
        telescope has periods)
        """
        coefficients = {name: getattr(self, name).tel[tel_id] for name in self.TERMS}
        periods = self._periods.get(tel_id, [])
        if len(periods) == 0:
            return coefficients
        if time is None:
            logger.warning("Pointing model of telescope %d: no time, coefficients of the parameters used", tel_id)
            return coefficients
        outside = np.ones(time.shape, dtype=bool)
        for start, stop, period_coefficients in periods:
            inside = np.ones(time.shape, dtype=bool)
            if start is not None:
                inside &= time >= start
            if stop is not None:
                inside &= time < stop
            outside &= ~inside
            for name in self.TERMS:
                coefficients[name] = np.where(inside, period_coefficients[name], coefficients[name])
        if np.any(outside):
            logger.warning(
                "Pointing model of telescope %d: %d positions outside of its periods, coefficients of the parameters used",
                tel_id, np.count_nonzero(outside),
            )
        return coefficients

    def corrections(self, tel_id, alt, az, time=None):
        alt = u.Quantity(alt, u.deg).to_value(u.rad)
        az = u.Quantity(az, u.deg).to_value(u.rad)
        c = self.coefficients(tel_id, time)
        with np.errstate(divide="ignore", invalid="ignore"):
            tan_alt = np.tan(alt)
            d_az = (
                -c["ia"] - c["an"] * np.sin(az) * tan_alt - c["aw"] * np.cos(az) * tan_alt
                - c["ca"] / np.cos(alt) - c["npae"] * tan_alt
                - c["aces"] * np.sin(az) + c["acec"] * np.cos(az)
            )
            d_alt = (
                c["ie"] - c["an"] * np.cos(az) + c["aw"] * np.sin(az) - c["tf"] * np.cos(alt) - c["tx"] / tan_alt
                + c["eces"] * np.sin(alt) + c["ecec"] * np.cos(alt)
            )
        return d_alt * u.arcsec, d_az * u.arcsec
