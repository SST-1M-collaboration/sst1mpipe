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
from ctapipe.containers import ArrayEventContainer, ReconstructedGeometryContainer
from ctapipe.reco.preprocessing import telescope_to_horizontal
from ctapipe.reco.reconstructor import ReconstructionProperty
from ctapipe.reco.stereo_combination import StereoMeanCombiner
from ctapipe.reco.telescope_event_handling import get_subarray_index
from ctapipe.reco.utils import add_defaults_and_meta

__all__ = ["StereoMARSCombiner", "mars_directions"]


def mars_directions(candidates, weights, groups, n_events):
    """
    MARS stereo direction of array events with two telescopes: the combination of the two possible
    positions of each telescope (signs of the disp) with the smallest angular distance, and the
    weighted mean of these positions.

    Only two telescopes: with one telescope, or more than two, the predicted signs are kept.

    Parameters
    ----------
    candidates: numpy.ndarray
        (n_telescope_events, 2, 3) unit vectors of the two positions, the predicted sign first
    weights, groups: numpy.ndarray
        (n_telescope_events,) weight and array event index of each telescope event
    n_events: int
        number of array events

    Returns
    -------
    alt, az, spread, distance: numpy.ndarray
        (n_events,) direction (rad), spread of the chosen positions, angular distance of the
        chosen positions of the two telescopes (rad, 0 otherwise); NaN without telescope event
    """
    # TODO: more than two telescopes, e.g. the combination of signs with the smallest sum of the
    #  pairwise distances (2 ** n_telescopes combinations)
    groups = np.asarray(groups, dtype=np.int64)  # ctapipe gives uint64 indices: not for bincount
    multiplicity = np.bincount(groups, minlength=n_events)
    chosen = candidates[:, 0].copy()
    distance = np.zeros(n_events)

    events = np.flatnonzero(multiplicity == 2)
    order = np.argsort(groups, kind="stable")
    start = np.cumsum(multiplicity) - multiplicity
    first, second = order[start[events]], order[start[events] + 1]
    # (n, 2, 2) distances of the positions of the first and second telescopes
    cos = np.einsum("nix,njx->nij", candidates[first], candidates[second])
    separation = np.arccos(np.clip(cos, -1, 1)).reshape(len(events), 4)
    best = np.argmin(separation, axis=1)
    chosen[first] = candidates[first, best // 2]
    chosen[second] = candidates[second, best % 2]
    distance[events] = separation[np.arange(len(events)), best]

    mean = np.zeros((n_events, 3))
    np.add.at(mean, groups, chosen * weights[:, np.newaxis])
    with np.errstate(invalid="ignore", divide="ignore"):
        mean /= np.bincount(groups, weights=weights, minlength=n_events)[:, np.newaxis]
        norm = np.linalg.norm(mean, axis=1)
        alt = np.arcsin(mean[:, 2] / norm)
        az = np.arctan2(mean[:, 1], mean[:, 0]) % (2 * np.pi)
        # https://en.wikipedia.org/wiki/Directional_statistics#Measures_of_location_and_spread
        spread = np.sqrt(-2 * np.log(np.clip(norm, 0, 1)))
    distance[multiplicity == 0] = np.nan
    return alt, az, spread, distance


def _candidates(disp, psi, fov_lon, fov_lat, pointing_alt, pointing_az):
    """(n, 2, 3) unit vectors of the positions at +disp (predicted) and -disp along the major axis"""
    disp, psi, fov_lon, fov_lat = (np.atleast_1d(x) for x in (disp, psi, fov_lon, fov_lat))
    vectors = []
    for sign in (1, -1):
        alt, az = telescope_to_horizontal(
            lon=fov_lon + sign * disp * np.cos(psi), lat=fov_lat + sign * disp * np.sin(psi),
            pointing_alt=pointing_alt, pointing_az=pointing_az,
        )
        alt, az = alt.to_value(u.rad), az.to_value(u.rad)
        vectors.append(np.stack([np.cos(alt) * np.cos(az), np.cos(alt) * np.sin(az), np.sin(alt)], axis=-1))
    return np.stack(vectors, axis=1)


class StereoMARSCombiner(StereoMeanCombiner):
    """
    Stereo combination of the MARS analysis, for two telescopes.

    Geometry of the disp reconstructor: each telescope gives two possible positions (its disp norm
    on either side of the centroid along the major axis); the combination of signs whose positions
    are the closest to each other is chosen, and the direction is their weighted mean (``weights``).
    ``goodness_of_fit``: angular distance of the chosen positions (deg). Energy and particle type:
    as the ``StereoMeanCombiner``.

    Not yet for more than two telescopes (TODO): events with one telescope, or more than two, keep
    the predicted signs (weighted mean of the predicted positions).

    Needs the telescope frame Hillas parameters and the telescope pointing.
    """

    def predict_table(self, mono_predictions: Table) -> Table:
        if self.property is not ReconstructionProperty.GEOMETRY:
            return super().predict_table(mono_predictions)

        obs_ids, event_ids, _, indices = get_subarray_index(mono_predictions)
        valid = np.asarray(mono_predictions[f"{self.prefix}_tel_is_valid"], dtype=bool)
        table = mono_predictions[valid]
        if "telescope_pointing_altitude" in table.colnames:
            pointing = table["telescope_pointing_altitude"].quantity, table["telescope_pointing_azimuth"].quantity
        else:
            pointing = table["subarray_pointing_lat"].quantity, table["subarray_pointing_lon"].quantity
        candidates = _candidates(
            table[f"{self.prefix}_tel_parameter"].quantity, table["hillas_psi"].quantity,
            table["hillas_fov_lon"].quantity, table["hillas_fov_lat"].quantity, *pointing,
        )
        weights = np.asarray(self._calculate_weights(table), dtype=float)
        alt, az, spread, distance = mars_directions(
            candidates, weights, indices[valid], len(obs_ids),
        )

        telescopes = [[] for _ in obs_ids]
        for index, tel_id in zip(indices[valid], table["tel_id"], strict=True):
            telescopes[index].append(tel_id)
        stereo = Table({
            "obs_id": obs_ids,
            "event_id": event_ids,
            f"{self.prefix}_alt": np.rad2deg(alt) * u.deg,
            f"{self.prefix}_az": np.rad2deg(az) * u.deg,
            f"{self.prefix}_ang_distance_uncert": np.rad2deg(spread) * u.deg,
            f"{self.prefix}_is_valid": np.isfinite(alt) & np.isfinite(az),
            f"{self.prefix}_goodness_of_fit": np.rad2deg(distance),
            f"{self.prefix}_telescopes": telescopes,
        })
        for name in ("obs_id", "event_id"):
            stereo[name].description = mono_predictions[name].description
        add_defaults_and_meta(stereo, ReconstructedGeometryContainer, self.prefix)
        return stereo

    def _combine_altaz(self, event: ArrayEventContainer):
        """Stereo direction of one array event (event loop)"""
        tel_ids, candidates, weights = [], [], []
        for tel_id, dl2 in event.dl2.tel.items():
            mono, disp = dl2.geometry.get(self.prefix), dl2.disp.get(self.prefix)
            if mono is None or disp is None or not mono.is_valid:
                continue
            dl1 = event.dl1.tel[tel_id].parameters
            pointing = event.monitoring.tel[tel_id].pointing
            candidates.append(_candidates(
                disp.parameter, dl1.hillas.psi, dl1.hillas.fov_lon, dl1.hillas.fov_lat,
                pointing.altitude, pointing.azimuth,
            )[0])
            weights.append(float(self._calculate_weights(dl1)))
            tel_ids.append(tel_id)

        if tel_ids:
            alt, az, spread, distance = (x[0] for x in mars_directions(
                np.array(candidates), np.array(weights), np.zeros(len(tel_ids)), 1,
            ))
        else:
            alt = az = spread = distance = np.nan
        event.dl2.stereo.geometry[self.prefix] = ReconstructedGeometryContainer(
            alt=np.rad2deg(alt) * u.deg,
            az=np.rad2deg(az) * u.deg,
            ang_distance_uncert=np.rad2deg(spread) * u.deg,
            goodness_of_fit=float(np.rad2deg(distance)),
            telescopes=tel_ids,
            is_valid=bool(tel_ids),
            prefix=self.prefix,
        )
