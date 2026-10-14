"""
Bending model of the SST-1M telescopes: offsets between the position of the drive (encoders) and
the direction pointed to in the sky, as a function of the altitude and azimuth. The drive adds
these offsets to the position of the target, so its logs give the position with the offsets
(`sst1mpipe.io.drive.read_drive_logs` reverts them with a `BendingModel`).
"""
import logging

import astropy.units as u
import numpy as np
from ctapipe.core import TelescopeComponent
from ctapipe.core.traits import Float, FloatTelescopeParameter, Int

__all__ = ["BendingModel"]

logger = logging.getLogger(__name__)


class BendingModel(TelescopeComponent):
    """
    Offsets (deg) added by the drive of a telescope to the altitude and azimuth of the direction
    pointed to in the sky: drive position = sky position + offsets(sky position), with

    - d_az = IA + CA / cos(alt) + NPAE tan(alt) + (AN sin(az) + AW cos(az)) tan(alt)
    - d_alt = IE + AN cos(az) - AW sin(az) + TF cos(alt) + TX / tan(alt)

    (the usual terms of the pointing models of alt-az mounts, with the signs written here).
    All the coefficients are 0 by default (no bending). A model of another form is a subclass
    defining `offsets`.

    The sky position of a drive position is found by iteration (`revert`), the offsets being
    given at the sky position.
    """

    IA = FloatTelescopeParameter(default_value=0.0, help="Azimuth index (deg)").tag(config=True)
    IE = FloatTelescopeParameter(default_value=0.0, help="Altitude index (deg)").tag(config=True)
    CA = FloatTelescopeParameter(
        default_value=0.0, help="Collimation: optical axis not perpendicular to the altitude axis (deg)"
    ).tag(config=True)
    NPAE = FloatTelescopeParameter(
        default_value=0.0, help="Non perpendicularity of the azimuth and altitude axes (deg)"
    ).tag(config=True)
    AN = FloatTelescopeParameter(default_value=0.0, help="Tilt of the azimuth axis to the north (deg)").tag(config=True)
    AW = FloatTelescopeParameter(default_value=0.0, help="Tilt of the azimuth axis to the west (deg)").tag(config=True)
    TF = FloatTelescopeParameter(default_value=0.0, help="Tube flexure, in cos(alt) (deg)").tag(config=True)
    TX = FloatTelescopeParameter(default_value=0.0, help="Tube flexure, in 1 / tan(alt) (deg)").tag(config=True)

    tolerance = Float(
        default_value=1e-6, help="Precision of the sky position found by revert (deg)"
    ).tag(config=True)
    max_iterations = Int(default_value=50, help="Maximum number of iterations of revert").tag(config=True)

    def offsets(self, tel_id, alt, az):
        """
        Offsets (d_alt, d_az) added by the drive of the telescope tel_id to the sky position
        alt, az (Quantity, angles)
        """
        alt = u.Quantity(alt, u.deg).to_value(u.rad)
        az = u.Quantity(az, u.deg).to_value(u.rad)
        c = {name: getattr(self, name).tel[tel_id] for name in ("IA", "IE", "CA", "NPAE", "AN", "AW", "TF", "TX")}
        with np.errstate(divide="ignore", invalid="ignore"):
            tan_alt = np.tan(alt)
            d_az = c["IA"] + c["CA"] / np.cos(alt) + c["NPAE"] * tan_alt + (c["AN"] * np.sin(az) + c["AW"] * np.cos(az)) * tan_alt
            d_alt = c["IE"] + c["AN"] * np.cos(az) - c["AW"] * np.sin(az) + c["TF"] * np.cos(alt) + c["TX"] / tan_alt
        return d_alt * u.deg, d_az * u.deg

    def apply(self, tel_id, alt, az):
        """Drive position (alt, az) of the telescope tel_id pointing to the sky position alt, az"""
        d_alt, d_az = self.offsets(tel_id, alt, az)
        return u.Quantity(alt, u.deg) + d_alt, u.Quantity(az, u.deg) + d_az

    def revert(self, tel_id, alt, az):
        """
        Sky position (alt, az) pointed to by the telescope tel_id at the drive position alt, az:
        inverse of `apply`, by fixed point iteration. The positions where the iteration does not
        converge (e.g. close to the zenith for large CA or NPAE) are NaN.
        """
        drive_alt = u.Quantity(alt, u.deg).to_value(u.deg)
        drive_az = u.Quantity(az, u.deg).to_value(u.deg)
        sky_alt, sky_az = drive_alt, drive_az
        for _ in range(self.max_iterations):
            d_alt, d_az = self.offsets(tel_id, sky_alt * u.deg, sky_az * u.deg)
            new_alt, new_az = drive_alt - d_alt.to_value(u.deg), drive_az - d_az.to_value(u.deg)
            change = np.maximum(np.abs(new_alt - sky_alt), np.abs(new_az - sky_az))
            sky_alt, sky_az = new_alt, new_az
            if np.all(change < self.tolerance):
                break
        else:
            failed = ~(change < self.tolerance)
            logger.warning(
                "Bending model of telescope %d: no convergence after %d iterations for %d positions (set to NaN)",
                tel_id, self.max_iterations, np.count_nonzero(failed),
            )
            sky_alt, sky_az = np.where(failed, np.nan, sky_alt), np.where(failed, np.nan, sky_az)
        return sky_alt * u.deg, sky_az * u.deg
