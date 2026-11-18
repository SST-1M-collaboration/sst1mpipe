"""
Stereo direction of the MARS analysis for the ctapipe disp reconstructor, selected in its configuration::

    DispReconstructor:
      stereo_combiner_cls: StereoMARSCombiner
      StereoMARSCombiner:
        weights: intensity

ctapipe finds the combiners among the imported subclasses of ``StereoCombiner``: ``sst1mpipe.ml`` must
be imported, as by the sst1mpipe tools.
"""
import astropy.units as u
import numpy as np
from astropy.table import Table
from ctapipe.containers import ArrayEventContainer
from ctapipe.reco.preprocessing import telescope_to_horizontal
from ctapipe.reco.reconstructor import ReconstructionProperty
from ctapipe.reco.stereo_combination import StereoMeanCombiner
from ctapipe.reco.telescope_event_handling import get_subarray_index

__all__ = ["StereoMARSCombiner", "mars_disp_signs"]


def mars_disp_signs(hillas, disps):
    """
    Signs of the disp of two telescopes (MARS): the crossing of the major axes of the two images is
    computed, and each telescope takes the sign of its disp whose position (centroid +/- disp along
    the major axis) is the closest to this crossing.

    The images are in the telescope frame: both telescopes are assumed to point to the same
    direction (same frame for both).

    Parameters
    ----------
    hillas: list of dict
        Hillas parameters of the telescopes: ``fov_lon``, ``fov_lat``, ``psi`` (quantities, or deg)
    disps: list
        disp of the telescopes (signed, quantity or deg)

    Returns
    -------
    numpy.ndarray
        corrected sign (+1 or -1) of the disp of each telescope; the predicted signs if the major
        axes are parallel (no crossing)

    Raises
    ------
    ValueError
        less than two telescopes
    NotImplementedError
        more than two telescopes (TODO)
    """
    if len(hillas) < 2 or len(disps) != len(hillas):
        raise ValueError(f"MARS disp signs: two telescopes needed, got {len(hillas)} Hillas and {len(disps)} disps")
    if len(hillas) > 2:
        raise NotImplementedError("MARS disp signs: more than two telescopes is not implemented yet")

    def deg(value):
        return u.Quantity(value, u.deg).to_value(u.deg)

    centroids = np.array([[deg(h["fov_lon"]), deg(h["fov_lat"])] for h in hillas])
    psi = np.deg2rad([deg(h["psi"]) for h in hillas])
    axes = np.stack([np.cos(psi), np.sin(psi)], axis=1)
    norms = np.abs([deg(d) for d in disps])
    predicted = np.where(np.array([deg(d) for d in disps]) < 0, -1, 1)

    # crossing: centroid_0 + t_0 axis_0 = centroid_1 + t_1 axis_1
    matrix = np.array([axes[0], -axes[1]]).T
    if abs(np.linalg.det(matrix)) < 1e-9:
        return predicted
    t = np.linalg.solve(matrix, centroids[1] - centroids[0])
    crossing = centroids[0] + t[0] * axes[0]

    signs = []
    for centroid, axis, norm in zip(centroids, axes, norms, strict=True):
        distance = {sign: np.linalg.norm(centroid + sign * norm * axis - crossing) for sign in (1, -1)}
        signs.append(min(distance, key=distance.get))
    return np.array(signs)


def _direction(fov_lon, fov_lat, psi, disp, pointing_alt, pointing_az):
    """alt, az of the position at disp from the centroid along the major axis"""
    return telescope_to_horizontal(
        lon=fov_lon + disp * np.cos(psi), lat=fov_lat + disp * np.sin(psi),
        pointing_alt=pointing_alt, pointing_az=pointing_az,
    )


class StereoMARSCombiner(StereoMeanCombiner):
    """
    Stereo combination of the MARS analysis: the ``StereoMeanCombiner`` with the disp signs of the
    telescopes corrected by `mars_disp_signs` (crossing of the major axes) for the events with two
    telescopes. The directions of the telescopes are recomputed with these signs, then combined by
    the ``StereoMeanCombiner`` (weighted mean, ``weights``). Events with one telescope keep the
    predicted sign. Energy and particle type: as the ``StereoMeanCombiner``.

    More than two telescopes: not implemented yet (TODO), raises ``NotImplementedError``.

    Needs the telescope frame Hillas parameters and the telescope pointing.
    """

    def predict_table(self, mono_predictions: Table) -> Table:
        if self.property is ReconstructionProperty.GEOMETRY:
            mono_predictions = self._correct_disp_signs(mono_predictions)
        return super().predict_table(mono_predictions)

    def _correct_disp_signs(self, mono_predictions):
        """Directions of the telescopes of the events with two valid telescopes, with the MARS disp signs"""
        table = mono_predictions.copy()
        prefix = f"{self.prefix}_tel"
        valid = np.asarray(table[f"{prefix}_is_valid"], dtype=bool)
        index = np.asarray(get_subarray_index(table)[3], dtype=np.int64)
        n_valid = np.bincount(index[valid], minlength=len(index))
        rows = np.flatnonzero(valid & (n_valid[index] >= 2))
        if len(rows) == 0:
            return table

        disp = table[f"{prefix}_parameter"].quantity[rows]
        for event in np.unique(index[rows]):
            event_rows = np.flatnonzero(index[rows] == event)
            hillas = [{name: table[f"hillas_{name}"].quantity[rows[i]] for name in ("fov_lon", "fov_lat", "psi")}
                      for i in event_rows]
            disp[event_rows] = mars_disp_signs(hillas, disp[event_rows]) * np.abs(disp[event_rows])

        if "telescope_pointing_altitude" in table.colnames:
            pointing = table["telescope_pointing_altitude"], table["telescope_pointing_azimuth"]
        else:
            pointing = table["subarray_pointing_lat"], table["subarray_pointing_lon"]
        alt, az = _direction(
            table["hillas_fov_lon"].quantity[rows], table["hillas_fov_lat"].quantity[rows],
            table["hillas_psi"].quantity[rows], disp, pointing[0].quantity[rows], pointing[1].quantity[rows],
        )
        table[f"{prefix}_alt"][rows] = alt.to_value(table[f"{prefix}_alt"].unit)
        table[f"{prefix}_az"][rows] = az.to_value(table[f"{prefix}_az"].unit)
        return table

    def _combine_altaz(self, event: ArrayEventContainer):
        """Stereo direction of one array event (event loop), as the table version"""
        tel_ids = [
            tel_id for tel_id, dl2 in event.dl2.tel.items()
            if self.prefix in dl2.geometry and self.prefix in dl2.disp and dl2.geometry[self.prefix].is_valid
        ]
        if len(tel_ids) < 2:
            super()._combine_altaz(event)
            return

        hillas = [event.dl1.tel[tel_id].parameters.hillas for tel_id in tel_ids]
        disps = [event.dl2.tel[tel_id].disp[self.prefix].parameter for tel_id in tel_ids]
        signs = mars_disp_signs([dict(fov_lon=h.fov_lon, fov_lat=h.fov_lat, psi=h.psi) for h in hillas], disps)
        geometries = [event.dl2.tel[tel_id].geometry[self.prefix] for tel_id in tel_ids]
        predicted = [(geometry.alt, geometry.az) for geometry in geometries]
        for tel_id, h, disp, sign, geometry in zip(tel_ids, hillas, disps, signs, geometries, strict=True):
            pointing = event.monitoring.tel[tel_id].pointing
            geometry.alt, geometry.az = _direction(h.fov_lon, h.fov_lat, h.psi, sign * abs(disp),
                                                   pointing.altitude, pointing.azimuth)
        try:
            super()._combine_altaz(event)
        finally:
            # the directions of the telescopes keep their predicted signs, as in the table version
            for geometry, (alt, az) in zip(geometries, predicted, strict=True):
                geometry.alt, geometry.az = alt, az
