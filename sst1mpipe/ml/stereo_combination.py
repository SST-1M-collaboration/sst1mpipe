"""
Stereo combination of the telescope predictions of the ctapipe reconstructors (``ctapipe.reco``)
for the SST-1M telescopes.

`StereoMARSCombiner` is selected in the configuration of a ctapipe reconstructor, e.g. of the
``DispReconstructor`` of ``sst1mpipe-train-disp-reconstructor``::

    DispReconstructor:
      stereo_combiner_cls: StereoMARSCombiner
      StereoMARSCombiner:
        weights: intensity

ctapipe finds the combiners among the imported subclasses of ``StereoCombiner``: ``sst1mpipe.ml``
must be imported, as by the sst1mpipe tools (``sst1mpipe-apply-models``, the training tools).
"""
import itertools

import astropy.units as u
import numpy as np
from astropy.table import Table
from ctapipe.compat import COPY_IF_NEEDED
from ctapipe.containers import ArrayEventContainer, ReconstructedGeometryContainer
from ctapipe.core.traits import Int
from ctapipe.reco.preprocessing import telescope_to_horizontal
from ctapipe.reco.reconstructor import ReconstructionProperty
from ctapipe.reco.stereo_combination import StereoMeanCombiner
from ctapipe.reco.telescope_event_handling import get_subarray_index
from ctapipe.reco.utils import add_defaults_and_meta

__all__ = ["StereoMARSCombiner", "mars_directions"]


def _unit_vectors(alt, az):
    """Cartesian unit vectors (..., 3) of horizontal directions (rad)"""
    cos_alt = np.cos(alt)
    return np.stack([cos_alt * np.cos(az), cos_alt * np.sin(az), np.sin(alt)], axis=-1)


def mars_directions(candidates, weights, groups, n_events=None, max_telescopes=10):
    """
    Stereo direction of array events with the MARS method: each telescope gives two possible
    source positions (the two signs of its disp); the combination of signs whose positions are the
    closest to each other (smallest sum of the pairwise angular distances) is chosen, and the
    direction is the weighted mean of its positions.

    Parameters
    ----------
    candidates: numpy.ndarray
        (n_telescope_events, 2, 3): unit vectors of the two possible positions of each telescope
        event (the predicted sign of the disp first)
    weights: numpy.ndarray
        (n_telescope_events,): weight of each telescope event
    groups: numpy.ndarray
        (n_telescope_events,): index of the array event of each telescope event
    n_events: int or None
        number of array events (at least the largest index of ``groups`` + 1)
    max_telescopes: int
        array events with more telescopes keep the predicted sign of each telescope (the number of
        combinations is 2 ** n_telescopes)

    Returns
    -------
    mean: numpy.ndarray
        (n_array_events, 3): weighted mean (not normalised) of the chosen positions
    distance: numpy.ndarray
        (n_array_events,): sum of the pairwise angular distances of the chosen positions (rad),
        0 for a single telescope
    chosen: numpy.ndarray
        (n_telescope_events, 3): unit vector of the chosen position of each telescope event
    """
    groups = np.asarray(groups, dtype=np.int64)  # e.g. uint64 indices of ctapipe: not for bincount
    if n_events is None:
        n_events = int(groups.max()) + 1 if len(groups) else 0
    chosen = candidates[:, 0].copy()
    distance = np.zeros(n_events)
    multiplicity = np.bincount(groups, minlength=n_events)
    # telescope events sorted by array event: the events of a multiplicity are a (n, m) index array
    order = np.argsort(groups, kind="stable")
    starts = np.concatenate([[0], np.cumsum(multiplicity)[:-1]])
    for m in np.unique(multiplicity):
        if m < 2 or m > max_telescopes:
            continue
        events = np.flatnonzero(multiplicity == m)
        rows = order[starts[events][:, np.newaxis] + np.arange(m)]  # (n, m) telescope events
        signs = np.array(list(itertools.product((0, 1), repeat=m)))  # (2 ** m, m)
        # positions of each combination: (n, 2 ** m, m, 3)
        positions = candidates[rows[:, np.newaxis, :], signs[np.newaxis, :, :]]
        pair_i, pair_j = np.triu_indices(m, k=1)
        cos = np.einsum("nckx,nckx->nck", positions[:, :, pair_i], positions[:, :, pair_j])
        total = np.arccos(np.clip(cos, -1, 1)).sum(axis=-1)  # (n, 2 ** m)
        best = np.argmin(total, axis=1)
        distance[events] = total[np.arange(len(events)), best]
        chosen[rows] = positions[np.arange(len(events)), best]
    mean = np.zeros((n_events, 3))
    np.add.at(mean, groups, chosen * weights[:, np.newaxis])
    return mean, distance, chosen


class StereoMARSCombiner(StereoMeanCombiner):
    """
    Stereo combination of the telescope predictions as in the MARS analysis (and in the stereo
    reconstruction of sst1mpipe, ``sst1mpipe.reco.stereo_reconstruction``).

    - geometry of the ``DispReconstructor``: each telescope gives two possible source positions,
      at the norm of its disp on either side of the image centroid along its major axis. The
      combination of signs whose positions are the closest to each other (smallest sum of the
      pairwise angular distances: the 4 combinations of 2 telescopes) is chosen, and the stereo
      direction is the weighted mean of its positions (``weights``). The sign predicted by the disp
      sign classifier is only used for the events with one telescope. The ``goodness_of_fit`` is the
      sum of the angular distances of the chosen positions (deg), the ``ang_distance_uncert`` their
      spread, as for the ``StereoMeanCombiner``.
    - energy and particle type: (weighted) mean, as the ``StereoMeanCombiner``.

    The telescope tables need the DL1 parameters in the telescope frame (``hillas_fov_lon``,
    ``hillas_fov_lat``, ``hillas_psi``) and the telescope pointing.
    """

    max_telescopes = Int(
        default_value=10,
        help=(
            "Array events with more telescopes keep the predicted disp sign of each telescope"
            " (2 ** n_telescopes combinations)"
        ),
    ).tag(config=True)

    def predict_table(self, mono_predictions: Table) -> Table:
        if self.property is not ReconstructionProperty.GEOMETRY:
            return super().predict_table(mono_predictions)

        prefix = f"{self.prefix}_tel"
        obs_ids, event_ids, _, tel_to_array_indices = get_subarray_index(mono_predictions)
        n_array_events = len(obs_ids)
        stereo_table = Table({"obs_id": obs_ids, "event_id": event_ids})
        for colname in ("obs_id", "event_id"):
            stereo_table[colname].description = mono_predictions[colname].description

        valid = np.asarray(mono_predictions[f"{prefix}_is_valid"], dtype=bool)
        # a single valid telescope event: the mean direction is a position, its spread 0
        alt = np.full(n_array_events, np.nan)
        az = np.full(n_array_events, np.nan)
        spread = np.full(n_array_events, np.nan)
        distance = np.full(n_array_events, np.nan)
        if np.count_nonzero(valid) > 0:
            table = mono_predictions[valid]
            weights = np.asarray(self._calculate_weights(table), dtype=float)
            candidates = self._candidates(table, prefix)
            groups = np.asarray(tel_to_array_indices[valid], dtype=np.int64)
            mean, sum_distance, _ = mars_directions(candidates, weights, groups, n_array_events, self.max_telescopes)
            weight_sum = np.bincount(groups, weights=weights, minlength=n_array_events)
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = mean / weight_sum[:, np.newaxis]
                norm = np.linalg.norm(mean, axis=1)
                alt = np.arcsin(mean[:, 2] / norm)
                az = np.arctan2(mean[:, 1], mean[:, 0]) % (2 * np.pi)
                # https://en.wikipedia.org/wiki/Directional_statistics#Measures_of_location_and_spread
                spread = np.sqrt(-2 * np.log(np.clip(norm, 0, 1)))
            has_valid = np.bincount(groups, minlength=n_array_events) > 0
            distance = np.where(has_valid, sum_distance, np.nan)

        stereo_table[f"{self.prefix}_alt"] = u.Quantity(np.rad2deg(alt), u.deg, copy=COPY_IF_NEEDED)
        stereo_table[f"{self.prefix}_az"] = u.Quantity(np.rad2deg(az), u.deg, copy=COPY_IF_NEEDED)
        stereo_table[f"{self.prefix}_ang_distance_uncert"] = u.Quantity(
            np.rad2deg(spread), u.deg, copy=COPY_IF_NEEDED
        )
        stereo_table[f"{self.prefix}_is_valid"] = np.isfinite(alt) & np.isfinite(az)
        stereo_table[f"{self.prefix}_goodness_of_fit"] = np.rad2deg(distance)

        tel_ids = [[] for _ in range(n_array_events)]
        for index, tel_id in zip(tel_to_array_indices[valid], mono_predictions["tel_id"][valid], strict=True):
            tel_ids[index].append(tel_id)
        stereo_table[f"{self.prefix}_telescopes"] = tel_ids
        add_defaults_and_meta(stereo_table, ReconstructedGeometryContainer, self.prefix)
        return stereo_table

    def _candidates(self, table, prefix):
        """
        Unit vectors (n, 2, 3) of the two possible source positions of the telescope events: the
        predicted one (signed disp) and the other one (opposite sign)
        """
        disp = table[f"{self.prefix}_tel_parameter"].quantity.to_value(u.deg)
        psi = table["hillas_psi"].quantity.to_value(u.rad)
        lon = table["hillas_fov_lon"].quantity.to_value(u.deg)
        lat = table["hillas_fov_lat"].quantity.to_value(u.deg)
        if "telescope_pointing_altitude" in table.colnames:
            pointing_alt = table["telescope_pointing_altitude"].quantity
            pointing_az = table["telescope_pointing_azimuth"].quantity
        else:
            pointing_alt = table["subarray_pointing_lat"].quantity
            pointing_az = table["subarray_pointing_lon"].quantity
        vectors = []
        for sign in (1, -1):
            alt, az = telescope_to_horizontal(
                lon=(lon + sign * disp * np.cos(psi)) * u.deg,
                lat=(lat + sign * disp * np.sin(psi)) * u.deg,
                pointing_alt=pointing_alt,
                pointing_az=pointing_az,
            )
            vectors.append(_unit_vectors(alt.to_value(u.rad), az.to_value(u.rad)))
        return np.stack(vectors, axis=1)

    def _combine_altaz(self, event: ArrayEventContainer):
        """Stereo direction of one array event (event loop, e.g. sst1mpipe-process)"""
        tel_ids, candidates, weights = [], [], []
        for tel_id, dl2 in event.dl2.tel.items():
            mono = dl2.geometry.get(self.prefix)
            disp = dl2.disp.get(self.prefix)
            if mono is None or not mono.is_valid or disp is None:
                continue
            dl1 = event.dl1.tel[tel_id].parameters
            pointing = event.monitoring.tel[tel_id].pointing
            table = Table({
                f"{self.prefix}_tel_parameter": u.Quantity([disp.parameter]),
                "hillas_psi": u.Quantity([dl1.hillas.psi]),
                "hillas_fov_lon": u.Quantity([dl1.hillas.fov_lon]),
                "hillas_fov_lat": u.Quantity([dl1.hillas.fov_lat]),
                "telescope_pointing_altitude": u.Quantity([pointing.altitude]),
                "telescope_pointing_azimuth": u.Quantity([pointing.azimuth]),
            })
            candidates.append(self._candidates(table, f"{self.prefix}_tel")[0])
            weights.append(float(self._calculate_weights(dl1)) if dl1 else 1.0)
            tel_ids.append(tel_id)

        if len(tel_ids) > 0:
            mean, distance, _ = mars_directions(
                np.array(candidates), np.array(weights), np.zeros(len(tel_ids), dtype=int), 1, self.max_telescopes,
            )
            mean = mean[0] / np.sum(weights)
            norm = np.linalg.norm(mean)
            alt = u.Quantity(np.rad2deg(np.arcsin(mean[2] / norm)), u.deg)
            az = u.Quantity(np.rad2deg(np.arctan2(mean[1], mean[0]) % (2 * np.pi)), u.deg)
            uncert = u.Quantity(np.rad2deg(np.sqrt(-2 * np.log(min(norm, 1.0)))), u.deg)
            goodness = float(np.rad2deg(distance[0]))
            valid = True
        else:
            alt = az = uncert = u.Quantity(np.nan, u.deg)
            goodness = np.nan
            valid = False

        event.dl2.stereo.geometry[self.prefix] = ReconstructedGeometryContainer(
            alt=alt,
            az=az,
            ang_distance_uncert=uncert,
            goodness_of_fit=goodness,
            telescopes=tel_ids,
            is_valid=valid,
            prefix=self.prefix,
        )
