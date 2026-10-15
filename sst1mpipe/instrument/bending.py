"""
Bending models of the SST-1M telescopes: pointing errors between the position commanded to the
drive (direction in the sky) and the position of the drive (encoders), as a function of the
altitude and azimuth. The drive adds these errors to the position of the target, so its logs
(`sst1mpipe.io.drive.read_drive_logs`) give the drive position, and `BendingModel.revert` gives
the direction pointed to in the sky.
"""
import logging
from abc import abstractmethod

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import Float, FloatTelescopeParameter, Int

__all__ = ["BendingModel", "MeeksBendingModel"]

logger = logging.getLogger(__name__)


class BendingModel(TelescopeComponent):
    """
    Base class of the bending models: pointing errors (d_alt, d_az) of a telescope, drive
    position = sky position + pointing errors(sky position). A child class defines
    `pointing_errors`.
    """

    tolerance = Float(
        default_value=1e-6, help="Precision of the sky position found by revert (deg)"
    ).tag(config=True)
    max_iterations = Int(default_value=50, help="Maximum number of iterations of revert").tag(config=True)

    @abstractmethod
    def pointing_errors(self, tel_id, alt, az):
        """
        Pointing errors (d_alt, d_az) of the telescope tel_id commanded to the sky position
        alt, az (Quantity, angles)
        """

    def apply(self, tel_id, pointing):
        """
        Drive position (AltAz SkyCoord) of the telescope tel_id commanded to the sky position
        ``pointing`` (AltAz SkyCoord). Altitudes beyond +-90 deg are NaN.
        """
        d_alt, d_az = self.pointing_errors(tel_id, pointing.alt, pointing.az)
        alt = (pointing.alt + d_alt).to_value(u.deg)
        alt = np.where(np.abs(alt) <= 90, alt, np.nan)
        return SkyCoord(alt=alt * u.deg, az=pointing.az + d_az, frame=pointing.frame)

    def revert(self, tel_id, pointing):
        """
        Sky position (AltAz SkyCoord) pointed to by the telescope tel_id at the drive position
        ``pointing`` (AltAz SkyCoord): inverse of `apply`, by fixed point iteration. The
        positions where the iteration does not converge (e.g. close to the zenith) are NaN.
        """
        drive_alt, drive_az = pointing.alt.to_value(u.deg), pointing.az.to_value(u.deg)
        alt, az = drive_alt, drive_az
        for _ in range(self.max_iterations):
            d_alt, d_az = self.pointing_errors(tel_id, alt * u.deg, az * u.deg)
            new_alt, new_az = drive_alt - d_alt.to_value(u.deg), drive_az - d_az.to_value(u.deg)
            change = np.maximum(np.abs(new_alt - alt), np.abs(new_az - az))
            alt, az = new_alt, new_az
            if np.all(change < self.tolerance):
                break
        else:
            failed = ~(change < self.tolerance)
            logger.warning(
                "Bending model of telescope %d: no convergence after %d iterations for %d positions (set to NaN)",
                tel_id, self.max_iterations, np.count_nonzero(failed),
            )
            alt, az = np.where(failed, np.nan, alt), np.where(failed, np.nan, az)
        alt = np.where(np.abs(alt) <= 90, alt, np.nan)
        return SkyCoord(alt=alt * u.deg, az=az * u.deg, frame=pointing.frame)


class MeeksBendingModel(BendingModel):
    """
    Bending model of alt-az mounts of Meeks, Ball and Hull [1]_: pointing errors (eq. 1:
    observed - command, the observed being the drive position) at the commanded position
    (elevation h, azimuth A) from the tilt of the azimuth axis (eq. 3), the skew of the elevation
    axis (eq. 5), the collimation error (eq. 6), the gravitational deflection (eq. 7) and the
    biases of the encoders (eqs. 8 and 9):

    - d_az = -phi cos(kappa_A - A) tan(h) + epsilon tan(h) - delta / cos(h) + K_A
    - d_alt = -phi sin(kappa_A - A) + beta h + K_h

    with h in deg in the gravitational deflection beta h + K_h. All the parameters are 0 by
    default (no pointing error).

    References
    ----------
    .. [1] M. L. Meeks, J. A. Ball and A. B. Hull, "The pointing calibration of the Haystack
       antenna", IEEE Transactions on Antennas and Propagation, vol. AP-16, no. 6, pp. 746-751
       (1968), https://doi.org/10.1109/TAP.1968.1139287
    """

    azimuth_axis_tilt = FloatTelescopeParameter(
        default_value=0.0, help="phi: amount of the tilt of the azimuth axis (deg)"
    ).tag(config=True)
    level_axis_azimuth = FloatTelescopeParameter(
        default_value=0.0, help="kappa_A: azimuth of the level axis, the azimuth axis is tilted toward kappa_A + 90 deg (deg)"
    ).tag(config=True)
    elevation_axis_skew = FloatTelescopeParameter(
        default_value=0.0, help="epsilon: skew of the elevation axis (deg)"
    ).tag(config=True)
    collimation = FloatTelescopeParameter(default_value=0.0, help="delta: collimation error (deg)").tag(config=True)
    elevation_coefficient = FloatTelescopeParameter(
        default_value=0.0, help="beta: gravitational deflection per degree of elevation"
    ).tag(config=True)
    azimuth_bias = FloatTelescopeParameter(default_value=0.0, help="K_A: azimuth bias (deg)").tag(config=True)
    elevation_bias = FloatTelescopeParameter(default_value=0.0, help="K_h: elevation bias (deg)").tag(config=True)

    def pointing_errors(self, tel_id, alt, az):
        h = u.Quantity(alt, u.deg).to_value(u.deg)
        az = u.Quantity(az, u.deg).to_value(u.rad)
        phi = self.azimuth_axis_tilt.tel[tel_id]
        kappa = np.deg2rad(self.level_axis_azimuth.tel[tel_id])
        with np.errstate(divide="ignore", invalid="ignore"):
            tan_h = np.tan(np.deg2rad(h))
            d_az = (
                -phi * np.cos(kappa - az) * tan_h
                + self.elevation_axis_skew.tel[tel_id] * tan_h
                - self.collimation.tel[tel_id] / np.cos(np.deg2rad(h))
                + self.azimuth_bias.tel[tel_id]
            )
        d_alt = -phi * np.sin(kappa - az) + self.elevation_coefficient.tel[tel_id] * h + self.elevation_bias.tel[tel_id]
        return d_alt * u.deg, d_az * u.deg
