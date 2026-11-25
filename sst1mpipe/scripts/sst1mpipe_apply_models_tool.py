"""
``sst1mpipe-apply-models``: ``ctapipe-apply-models`` (same options) for the SST-1M files. The output
is the DL2 of ctapipe.

    sst1mpipe-apply-models -i run.dl1.h5 -o run.dl2.h5 \\
        -r models/nsb100/energy_regressor.pkl -r models/nsb100/particle_classifier.pkl \\
        -r models/nsb100/disp_reconstructor.pkl

The models use the pointing of the telescopes (``telescope_pointing_altitude``,
``telescope_pointing_azimuth``: e.g. the feature ``cos_zenith``, and the direction of the disp). Files
without it (data processed without monitoring data) use the pointing of the target, see
`add_telescope_pointing`.
"""
from unittest import mock

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time
from ctapipe.containers import CoordinateFrameType
from ctapipe.core.traits import Int, IntTelescopeParameter, List
from ctapipe.reco import Reconstructor
from ctapipe.tools.apply_models import ApplyModels as CtapipeApplyModels

import sst1mpipe.ml  # noqa: F401 (stereo combiners of sst1mpipe)

#: Reconstructor.read of ctapipe (replaced during ApplyModels.setup)
_read_reconstructor = Reconstructor.read

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
    pointing of the telescope (warning), see `add_telescope_pointing`, and optionally only for some
    telescopes (``telescopes``: e.g. the mono analysis of one telescope).

    The models are trained per telescope type of the simulations: each telescope of the input file
    uses the models of the simulated telescope ``model_telescopes`` (data telescopes 21, 22: the
    simulated telescopes 1, 2).
    """

    name = "sst1mpipe-apply-models"

    telescopes = List(
        Int(),
        default_value=[],
        help=(
            "Telescopes to which the models are applied (all if empty): e.g. [1] for the analysis of"
            " telescope 1 alone, the DL2 of the subarray is then its mono reconstruction. The other telescopes"
            " get invalid predictions (no Hillas parameters: the models must use one of them)"
        ),
    ).tag(config=True)

    model_telescopes = IntTelescopeParameter(
        default_value=[("type", "*", 0), ("id", 21, 1), ("id", 22, 2)],
        help="Telescope of the models applied to each telescope of the input file (0: the same telescope id)",
    ).tag(config=True)

    def setup(self):
        with mock.patch.object(Reconstructor, "read", self._read_reconstructor):
            super().setup()

    def _read_reconstructor(self, path, parent=None, subarray=None, **kwargs):
        """The models of the telescopes of the simulations, for the telescopes of the input file"""
        reconstructor = _read_reconstructor(path, parent=parent, **kwargs)
        models, model_subarray = reconstructor._models, reconstructor.subarray
        self.model_telescopes.attach_subarray(subarray)
        reconstructor._models = {}
        for tel_id, tel in subarray.tel.items():
            model_tel_id = self.model_telescopes.tel[tel_id] or tel_id
            model_tel = model_subarray.tel.get(model_tel_id)
            if model_tel is not None and model_tel in models:
                if reconstructor._models.get(tel, models[model_tel]) is not models[model_tel]:
                    raise ValueError(f"Telescope {tel_id}: same telescope type as another telescope, other models")
                reconstructor._models[tel] = models[model_tel]
                self.log.info("%s: models of telescope %d for telescope %d", type(reconstructor).__name__,
                              model_tel_id, tel_id)
        reconstructor.subarray = subarray
        return reconstructor

    def _apply(self, reconstructor, tel_tables, start, stop):
        if self.telescopes:
            # the other telescopes: no Hillas parameters, invalid predictions (each subarray event keeps
            # its row in the DL2 tables, invalid without the selected telescopes)
            for tel_id, table in tel_tables.items():
                if tel_id not in self.telescopes:
                    for name in table.colnames:
                        if name.startswith("hillas_"):
                            nan = np.full(len(table), np.nan)
                            unit = getattr(table[name], "unit", None)
                            table[name] = nan if unit is None else u.Quantity(nan, unit)
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
