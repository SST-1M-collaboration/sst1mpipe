"""
Stereo direction of the ctapipe disp reconstructor from the disp norms of the telescopes and a choice
of their disp signs (`StereoDispSignCombiner`):

- `StereoDispSignClosestCombiner`: the SST-1M Crab paper (arXiv:2506.01733, Sec. 4.3.2, the MAGIC
  approach, as sst1mpipe < 0.9), the closest combination of the signs;
- `StereoDispSignCrossingCombiner`: the former ``StereoMARSCombiner`` of sst1mpipe, the signs closest
  to the crossing of the major axes.

Selected in the configuration of the disp reconstructor::

    DispReconstructor:
      stereo_combiner_cls: StereoDispSignClosestCombiner
      StereoDispSignClosestCombiner:
        weights: intensity
        min_telescopes: 2

or when the models are applied (``sst1mpipe-apply-models --ApplyModels.disp_stereo_combiner=<name>``).

ctapipe finds the combiners among the imported subclasses of ``StereoCombiner``: ``sst1mpipe.ml`` must
be imported, as by the sst1mpipe tools.
"""
from abc import abstractmethod
from itertools import product

import astropy.units as u
import numpy as np
from astropy.table import Table
from ctapipe.containers import ArrayEventContainer
from ctapipe.core.traits import CaselessStrEnum, Int
from ctapipe.reco.preprocessing import telescope_to_horizontal
from ctapipe.reco.reconstructor import ReconstructionProperty
from ctapipe.reco.stereo_combination import StereoMeanCombiner
from ctapipe.reco.telescope_event_handling import get_subarray_index

__all__ = [
    "StereoDispSignCombiner",
    "StereoDispSignClosestCombiner",
    "StereoDispSignCrossingCombiner",
    "disp_sign_combination",
    "axes_crossing_disp_signs",
]


def _unit_vectors(alt, az):
    alt, az = u.Quantity(alt, u.rad).value, u.Quantity(az, u.rad).value
    return np.stack([np.cos(alt) * np.cos(az), np.cos(alt) * np.sin(az), np.sin(alt)], axis=-1)


def disp_sign_combination(alt_plus, az_plus, alt_minus, az_minus):
    """
    Disp signs of the telescopes of events (SST-1M Crab paper, MAGIC): each telescope gives two
    possible source positions, at its disp norm on either side of the image centroid along the major
    axis (signs +1 and -1). Of all the combinations of the signs of the telescopes (4 for two
    telescopes), the one whose positions are the closest is chosen: for two telescopes the smallest
    angular distance between the two positions, for more telescopes the smallest sum of the angular
    distances of all the pairs.

    Parameters
    ----------
    alt_plus, az_plus, alt_minus, az_minus: Quantity, shape (n_events, n_telescopes)
        positions of the telescopes of the events with the signs +1 and -1

    Returns
    -------
    numpy.ndarray, shape (n_events, n_telescopes)
        sign (+1 or -1) of the disp of each telescope

    Raises
    ------
    ValueError
        less than two telescopes
    """
    plus, minus = _unit_vectors(alt_plus, az_plus), _unit_vectors(alt_minus, az_minus)
    if plus.ndim == 2:  # one event
        return disp_sign_combination(*(np.atleast_2d(x) for x in (alt_plus, az_plus, alt_minus, az_minus)))[0]
    n_events, n_tel = plus.shape[:2]
    if n_tel < 2:
        raise ValueError(f"disp sign combination: two telescopes needed, got {n_tel}")
    combinations = np.array(list(product((1, -1), repeat=n_tel)))  # (n_combinations, n_tel)
    pairs = [(i, j) for i in range(n_tel) for j in range(i + 1, n_tel)]
    distance = np.zeros((len(combinations), n_events))
    for c, signs in enumerate(combinations):
        vectors = np.where((signs > 0)[None, :, None], plus, minus)
        for i, j in pairs:
            cos = np.clip(np.sum(vectors[:, i] * vectors[:, j], axis=-1), -1, 1)
            distance[c] += np.arccos(cos)
    best = np.argmin(np.nan_to_num(distance, nan=np.inf), axis=0)
    return combinations[best]


def axes_crossing_disp_signs(hillas, disps):
    """
    Signs of the disp of two telescopes by the crossing of the major axes (`StereoDispSignCrossingCombiner`,
    the former ``StereoMARSCombiner`` of sst1mpipe; not the method of the SST-1M Crab paper): the crossing of the major axes of the two images is
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
        raise ValueError(f"axes crossing disp signs: two telescopes needed, got {len(hillas)} Hillas and {len(disps)} disps")
    if len(hillas) > 2:
        raise NotImplementedError("axes crossing disp signs: more than two telescopes is not implemented yet")

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


def _positions(fov_lon, fov_lat, psi, disp_norm, pointing_alt, pointing_az):
    """alt, az of the two positions at disp_norm from the centroid along the major axis (signs +1, -1)"""
    return [
        telescope_to_horizontal(
            lon=fov_lon + sign * disp_norm * np.cos(psi), lat=fov_lat + sign * disp_norm * np.sin(psi),
            pointing_alt=pointing_alt, pointing_az=pointing_az,
        )
        for sign in (1, -1)
    ]


class StereoDispSignCombiner(StereoMeanCombiner):
    """
    Base class of the stereo combiners of the direction which choose the disp signs of the telescopes
    from their disp norms (not their sign classifier) for the events with at least two valid telescopes
    (`_signs`): each telescope gives two possible source positions, at its disp norm on either side of
    the image centroid along the major axis. The direction is the weighted mean of the positions of the
    chosen signs (as the ``StereoMeanCombiner``, ``weights``: intensity by default, as the paper).

    Stereo events are the events with at least ``min_telescopes`` valid telescopes; the other ones get
    an invalid direction. The events with one valid telescope (``min_telescopes`` 1) keep the sign of
    the classifier.

    Needs the telescope frame Hillas parameters and the telescope pointing. Energy and particle type:
    as the ``StereoMeanCombiner``.
    """

    weights = CaselessStrEnum(
        ["none", "intensity", "aspect-weighted-intensity", "containment-weighted-intensity"],
        default_value="intensity",
        help="Weights of the mean of the directions of the telescopes (see StereoMeanCombiner); paper: intensity",
    ).tag(config=True)

    min_telescopes = Int(
        default_value=2, min=1,
        help=(
            "Minimum number of valid telescopes of a stereo event (SST-1M Crab paper: 2, both telescopes pass"
            " the quality cuts); the events with fewer telescopes get an invalid direction"
        ),
    ).tag(config=True)

    @abstractmethod
    def _signs(self, hillas, disps, plus, minus):
        """
        Disp signs (+1, -1) of the telescopes of events with the same number of valid telescopes n.

        Parameters
        ----------
        hillas: dict
            ``fov_lon``, ``fov_lat``, ``psi`` (Quantity, shape (n_events, n))
        disps: Quantity, shape (n_events, n)
            predicted (signed) disp
        plus, minus: tuple of Quantity
            (alt, az) of the positions with the signs +1 and -1, shape (n_events, n)

        Returns
        -------
        numpy.ndarray, shape (n_events, n)
        """

    def predict_table(self, mono_predictions: Table) -> Table:
        if self.property is not ReconstructionProperty.GEOMETRY:
            return super().predict_table(mono_predictions)
        table = self._correct_disp_signs(mono_predictions)
        stereo = super().predict_table(table)
        valid = np.asarray(table[f"{self.prefix}_tel_is_valid"], dtype=bool)
        index = np.asarray(get_subarray_index(table)[3], dtype=np.int64)
        n_valid = np.bincount(index[valid], minlength=len(stereo))
        too_few = n_valid < self.min_telescopes
        for name in (f"{self.prefix}_alt", f"{self.prefix}_az", f"{self.prefix}_ang_distance_uncert"):
            stereo[name][too_few] = np.nan
        stereo[f"{self.prefix}_is_valid"] = np.asarray(stereo[f"{self.prefix}_is_valid"], dtype=bool) & ~too_few
        return stereo

    def _correct_disp_signs(self, mono_predictions):
        """Directions of the telescopes of the events with at least two valid telescopes, with the chosen signs"""
        table = mono_predictions.copy()
        prefix = f"{self.prefix}_tel"
        valid = np.asarray(table[f"{prefix}_is_valid"], dtype=bool)
        index = np.asarray(get_subarray_index(table)[3], dtype=np.int64)
        n_valid = np.bincount(index[valid], minlength=index.max() + 1 if len(index) else 0)
        rows = np.flatnonzero(valid & (n_valid[index] >= 2))
        if len(rows) == 0:
            return table

        if "telescope_pointing_altitude" in table.colnames:
            pointing = table["telescope_pointing_altitude"].quantity, table["telescope_pointing_azimuth"].quantity
        else:
            pointing = table["subarray_pointing_lat"].quantity, table["subarray_pointing_lon"].quantity
        hillas = {name: table[f"hillas_{name}"].quantity[rows] for name in ("fov_lon", "fov_lat", "psi")}
        disp = table[f"{prefix}_parameter"].quantity[rows]
        plus, minus = _positions(hillas["fov_lon"], hillas["fov_lat"], hillas["psi"], np.abs(disp),
                                 pointing[0][rows], pointing[1][rows])
        signs = np.ones(len(rows), dtype=int)
        # events grouped by their number of valid telescopes, the rows of an event together
        multiplicity = n_valid[index[rows]]
        for n_tel in np.unique(multiplicity):
            selected = np.flatnonzero(multiplicity == n_tel)
            order = selected[np.argsort(index[rows][selected], kind="stable")]

            def grouped(values, order=order, n_tel=n_tel):
                return values[order].reshape(-1, n_tel)

            signs[order] = self._signs(
                {name: grouped(value) for name, value in hillas.items()}, grouped(disp),
                (grouped(plus[0]), grouped(plus[1])), (grouped(minus[0]), grouped(minus[1])),
            ).ravel()
        alt = np.where(signs > 0, plus[0].to_value(u.deg), minus[0].to_value(u.deg)) * u.deg
        az = np.where(signs > 0, plus[1].to_value(u.deg), minus[1].to_value(u.deg)) * u.deg
        table[f"{prefix}_alt"][rows] = alt.to_value(table[f"{prefix}_alt"].unit or u.deg)
        table[f"{prefix}_az"][rows] = az.to_value(table[f"{prefix}_az"].unit or u.deg)
        return table

    def _combine_altaz(self, event: ArrayEventContainer):
        """Stereo direction of one array event (event loop), as the table version"""
        tel_ids = [
            tel_id for tel_id, dl2 in event.dl2.tel.items()
            if self.prefix in dl2.geometry and self.prefix in dl2.disp and dl2.geometry[self.prefix].is_valid
        ]
        if len(tel_ids) < 2:
            super()._combine_altaz(event)
            self._invalidate_too_few(event, len(tel_ids))
            return

        parameters = [event.dl1.tel[tel_id].parameters.hillas for tel_id in tel_ids]
        def row(values):
            """quantities of the telescopes, shape (1, n_telescopes)"""
            return u.Quantity(values).reshape(1, -1)

        hillas = {name: row([getattr(h, name) for h in parameters]) for name in ("fov_lon", "fov_lat", "psi")}
        disp = row([event.dl2.tel[tel_id].disp[self.prefix].parameter for tel_id in tel_ids])
        pointing = [event.monitoring.tel[tel_id].pointing for tel_id in tel_ids]
        plus, minus = _positions(hillas["fov_lon"], hillas["fov_lat"], hillas["psi"], np.abs(disp),
                                 row([p.altitude for p in pointing]), row([p.azimuth for p in pointing]))
        signs = self._signs(hillas, disp, plus, minus)[0]
        geometries = [event.dl2.tel[tel_id].geometry[self.prefix] for tel_id in tel_ids]
        predicted = [(geometry.alt, geometry.az) for geometry in geometries]
        for i, (geometry, sign) in enumerate(zip(geometries, signs, strict=True)):
            position = plus if sign > 0 else minus
            geometry.alt, geometry.az = position[0][0, i], position[1][0, i]
        try:
            super()._combine_altaz(event)
        finally:
            # the directions of the telescopes keep their predicted signs, as in the table version
            for geometry, (alt, az) in zip(geometries, predicted, strict=True):
                geometry.alt, geometry.az = alt, az
        self._invalidate_too_few(event, len(tel_ids))

    def _invalidate_too_few(self, event, n_telescopes):
        """invalid stereo direction for less than ``min_telescopes`` valid telescopes"""
        if n_telescopes < self.min_telescopes:
            geometry = event.dl2.stereo.geometry[self.prefix]
            geometry.alt, geometry.az = np.nan * u.deg, np.nan * u.deg
            geometry.ang_distance_uncert = np.nan * u.deg
            geometry.is_valid = False


class StereoDispSignClosestCombiner(StereoDispSignCombiner):
    """
    Stereo direction of the SST-1M Crab paper (arXiv:2506.01733, Sec. 4.3.2; the MAGIC approach,
    ``stereo_reconstruction`` of sst1mpipe < 0.9): of the combinations of the disp signs of the
    telescopes, the one whose positions are the closest (`disp_sign_combination`); intensity weighted
    mean; stereo events: both telescopes valid (``min_telescopes`` 2). See `StereoDispSignCombiner`.
    """

    def _signs(self, hillas, disps, plus, minus):
        return disp_sign_combination(plus[0], plus[1], minus[0], minus[1])


class StereoDispSignCrossingCombiner(StereoDispSignCombiner):
    """
    The former ``StereoMARSCombiner`` of sst1mpipe: each telescope takes the disp sign whose position is
    the closest to the crossing of the major axes of the two telescopes (`axes_crossing_disp_signs`,
    two telescopes only: ``NotImplementedError`` for more); intensity weighted mean. By default the
    events with one valid telescope are kept with the sign of its classifier (``min_telescopes`` 1, as
    the ``StereoMARSCombiner``). Not the method of the SST-1M Crab paper (`StereoDispSignClosestCombiner`).
    """

    min_telescopes = Int(
        default_value=1, min=1,
        help="Minimum number of valid telescopes of a stereo event, see StereoDispSignCombiner",
    ).tag(config=True)

    def _signs(self, hillas, disps, plus, minus):
        return np.array([
            axes_crossing_disp_signs(
                [{name: hillas[name][event, i] for name in hillas} for i in range(disps.shape[1])], disps[event],
            )
            for event in range(disps.shape[0])
        ])
