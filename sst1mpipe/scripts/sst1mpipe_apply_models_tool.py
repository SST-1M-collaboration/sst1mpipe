"""
``sst1mpipe-apply-models``: ``ctapipe-apply-models`` (same options) for the SST-1M files. The output
is the DL2 of ctapipe; the models can use the ``StereoMARSCombiner``.

    sst1mpipe-apply-models -i run.dl1.h5 -o run.dl2.h5 \\
        -r models/nsb100/energy_regressor.pkl -r models/nsb100/particle_classifier.pkl \\
        -r models/nsb100/disp_reconstructor.pkl

The models use the pointing of the telescopes (``telescope_pointing_altitude``,
``telescope_pointing_azimuth``: e.g. the feature ``cos_zenith``, and the direction of the disp). Files
without it (data processed without monitoring data) use the pointing of the target, see
`add_telescope_pointing`.
"""
import astropy.units as u
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time
from ctapipe.containers import CoordinateFrameType
from ctapipe.tools.apply_models import ApplyModels as CtapipeApplyModels

import sst1mpipe.ml  # noqa: F401 (StereoMARSCombiner)

__all__ = ["ApplyModels", "add_telescope_pointing"]


def add_telescope_pointing(table, location):
    """
    Pointing of the telescope (``telescope_pointing_altitude``, ``telescope_pointing_azimuth``) of
    the telescope events of ``table``, when it is not in the table: the pointing of the target (in
    ICRS: ``subarray_pointing_lon`` / ``lat``, transformed into alt/az at the time of the events and
    at ``location``; or in alt/az).

    Returns
    -------
    bool: whether the pointing of the target was added
    """
    if "telescope_pointing_altitude" in table.colnames or "subarray_pointing_frame" not in table.colnames:
        return False
    frame = CoordinateFrameType(int(table["subarray_pointing_frame"][0]))
    lon, lat = u.Quantity(table["subarray_pointing_lon"], u.deg), u.Quantity(table["subarray_pointing_lat"], u.deg)
    if frame is CoordinateFrameType.ICRS:
        time = table["time"] if isinstance(table["time"], Time) else Time(table["time"])
        altaz = SkyCoord(ra=lon, dec=lat, frame="icrs").transform_to(AltAz(obstime=time, location=location))
        altitude, azimuth = altaz.alt, altaz.az
    elif frame is CoordinateFrameType.ALTAZ:
        altitude, azimuth = lat, lon
    else:
        return False
    table["telescope_pointing_altitude"] = altitude.to(u.rad)
    table["telescope_pointing_azimuth"] = azimuth.to(u.rad)
    return True


class ApplyModels(CtapipeApplyModels):
    """
    ``ctapipe-apply-models`` with the pointing of the target for the telescope events without the
    pointing of the telescope (warning), see `add_telescope_pointing`.
    """

    name = "sst1mpipe-apply-models"

    def _apply(self, reconstructor, tel_tables, start, stop):
        subarray = self.loader.subarray
        for tel_id, table in tel_tables.items():
            location = subarray.tel_earth_locations.get(tel_id, subarray.reference_location)
            if add_telescope_pointing(table, location) and not getattr(self, "_warned_pointing", False):
                self.log.warning(
                    "No pointing of the telescopes in %s: pointing of the target used (alt/az at the time"
                    " of the events)", self.input_url,
                )
                self._warned_pointing = True
        super()._apply(reconstructor, tel_tables, start, stop)


def main():
    ApplyModels().run()
