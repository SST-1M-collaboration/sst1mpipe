"""
Bending models of the SST-1M telescopes: corrections between the position of the drive
(encoders) and the direction pointed to in the sky, as a function of the altitude and azimuth.
The drive applies its bending model to the position of the target, so its logs
(`sst1mpipe.io.drive.read_drive_logs`) give the drive position, and `BendingModel.revert` gives
the direction pointed to in the sky.
"""
import logging
from abc import abstractmethod

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import Float, FloatTelescopeParameter

__all__ = ["BendingModel", "TPointBendingModel"]

logger = logging.getLogger(__name__)


class BendingModel(TelescopeComponent):
    """
    Base class of the bending models: corrections (d_alt, d_az) of the drive position of a
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
    def corrections(self, tel_id, alt, az):
        """
        Corrections (d_alt, d_az) of the drive position alt, az (Quantity, angles) of the
        telescope tel_id
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
                "Bending model of telescope %d: %d positions inverted with an error above %.1f arcsec (largest %.1f arcsec)",
                tel_id, np.count_nonzero(above), self.tolerance, np.max(distance[above]),
            )
        return drive

    def _add_corrections(self, tel_id, pointing, sign):
        """``pointing`` + sign * corrections at ``pointing``"""
        d_alt, d_az = self.corrections(tel_id, pointing.alt, pointing.az)
        alt = (pointing.alt + sign * d_alt).to_value(u.deg)
        alt = np.where(np.abs(alt) <= 90, alt, np.nan)
        return SkyCoord(alt=alt * u.deg, az=pointing.az + sign * d_az, frame=pointing.frame)


class TPointBendingModel(BendingModel):
    """
    Bending model of alt-az mounts with the 8 terms of TPOINT [1]_ (names, signs and units), as
    implemented in pytelpoint [2]_ (``pytelpoint.transform.azel_model``), corrections of the drive
    (encoder) position alt, az:

    - d_az = -IA - AN sin(az) tan(alt) - AW cos(az) tan(alt) - CA / cos(alt) - NPAE tan(alt)
    - d_alt = IE - AN cos(az) + AW sin(az) - TF cos(alt) - TX / tan(alt)

    The coefficients are in arcsec, 0 by default (no correction), configured with the TPOINT names
    in lower case (``ia``, ``ie``, ``an``, ``aw``, ``ca``, ``npae``, ``tf``, ``tx``), e.g.
    ``{"TPointBendingModel": {"ia": [["type", "*", 0], ["id", 21, -8360], ["id", 22, -2070]]}}``.

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

    def corrections(self, tel_id, alt, az):
        alt = u.Quantity(alt, u.deg).to_value(u.rad)
        az = u.Quantity(az, u.deg).to_value(u.rad)
        c = {name: getattr(self, name).tel[tel_id] for name in ("ia", "ie", "an", "aw", "ca", "npae", "tf", "tx")}
        with np.errstate(divide="ignore", invalid="ignore"):
            tan_alt = np.tan(alt)
            d_az = (
                -c["ia"] - c["an"] * np.sin(az) * tan_alt - c["aw"] * np.cos(az) * tan_alt
                - c["ca"] / np.cos(alt) - c["npae"] * tan_alt
            )
            d_alt = c["ie"] - c["an"] * np.cos(az) + c["aw"] * np.sin(az) - c["tf"] * np.cos(alt) - c["tx"] / tan_alt
        return d_alt * u.arcsec, d_az * u.arcsec
