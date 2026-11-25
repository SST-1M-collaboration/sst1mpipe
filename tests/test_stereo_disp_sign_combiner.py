"""Stereo combiners of the direction from the disp norms: closest combination of the signs (paper), axes crossing"""
from itertools import product

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import angular_separation
from astropy.table import Table
from ctapipe.containers import (
    ArrayEventContainer,
    DispContainer,
    ImageParametersContainer,
    ReconstructedGeometryContainer,
)
from ctapipe.reco.preprocessing import telescope_to_horizontal
from ctapipe.reco.reconstructor import ReconstructionProperty
from ctapipe.reco.stereo_combination import StereoCombiner, StereoMeanCombiner

from sst1mpipe.ml import (
    StereoDispSignClosestCombiner,
    StereoDispSignCrossingCombiner,
    axes_crossing_disp_signs,
    disp_sign_combination,
)

GEOMETRY = ReconstructionProperty.GEOMETRY


def test_combiner_is_available():
    """the combiner can be chosen by name in the configuration of the reconstructors"""
    names = StereoCombiner.non_abstract_subclasses()
    assert {"StereoDispSignClosestCombiner", "StereoDispSignCrossingCombiner"} <= set(names)
    assert "StereoDispSignCombiner" not in names  # base class
    closest = StereoDispSignClosestCombiner(prefix="disp", property=GEOMETRY)
    assert closest.weights == "intensity" and closest.min_telescopes == 2
    crossing = StereoDispSignCrossingCombiner(prefix="disp", property=GEOMETRY)
    assert crossing.weights == "intensity" and crossing.min_telescopes == 1  # as the former StereoMARSCombiner


def test_disp_sign_combination():
    """the combination of the signs whose positions are the closest"""
    alt, az = u.deg, u.deg
    # tel 1: + at (70, 10), - at (71, 10); tel 2: + at (70.5, 10.5), - at (70.05, 10): closest (+, -)
    signs = disp_sign_combination([70, 70.5] * alt, [10, 10.5] * az, [71, 70.05] * alt, [10, 10] * az)
    np.testing.assert_array_equal(signs, [1, -1])
    # several events at once
    signs = disp_sign_combination([[70, 70.5], [71, 70.05]] * alt, [[10, 10.5], [10, 10]] * az,
                                  [[71, 70.05], [70, 70.5]] * alt, [[10, 10], [10, 10.5]] * az)
    np.testing.assert_array_equal(signs, [[1, -1], [-1, 1]])
    # three telescopes: smallest sum of the distances of the pairs
    signs = disp_sign_combination([70, 72, 70.1] * alt, [10, 10, 10] * az, [72, 70.05, 73] * alt, [10, 10, 10] * az)
    np.testing.assert_array_equal(signs, [1, -1, 1])
    with pytest.raises(ValueError, match="two telescopes"):
        disp_sign_combination([[70]] * alt, [[10]] * az, [[71]] * alt, [[10]] * az)


def direction(table, disp):
    psi = table["hillas_psi"].quantity
    alt, az = telescope_to_horizontal(
        lon=table["hillas_fov_lon"].quantity + disp * np.cos(psi),
        lat=table["hillas_fov_lat"].quantity + disp * np.sin(psi),
        pointing_alt=table["telescope_pointing_altitude"].quantity,
        pointing_az=table["telescope_pointing_azimuth"].quantity,
    )
    return alt.to(u.deg), az.to(u.deg)


def mono_table(n_events=60, n_tel=2, seed=0):
    """telescope events with the columns of the disp reconstructor (pointing 70 deg altitude)"""
    rng = np.random.default_rng(seed)
    rows = []
    for event_id in range(n_events):
        for tel_id in range(1, n_tel + 1):
            rows.append(dict(
                obs_id=1, event_id=event_id, tel_id=tel_id,
                hillas_fov_lon=rng.uniform(-1, 1), hillas_fov_lat=rng.uniform(-1, 1),
                hillas_psi=rng.uniform(-90, 90), hillas_intensity=rng.uniform(50, 5000),
                disp_tel_parameter=rng.uniform(-1.5, 1.5), disp_tel_is_valid=rng.uniform() > 0.2,
                telescope_pointing_altitude=70.0, telescope_pointing_azimuth=10.0,
            ))
    table = Table(rows=rows)
    for name in ("hillas_fov_lon", "hillas_fov_lat", "hillas_psi", "disp_tel_parameter",
                 "telescope_pointing_altitude", "telescope_pointing_azimuth"):
        table[name].unit = u.deg
    table["disp_tel_alt"], table["disp_tel_az"] = direction(table, table["disp_tel_parameter"].quantity)
    return table


def brute_force_signs(table, rows):
    """signs of the closest positions, by computing the distances of all the combinations"""
    norm = np.abs(table["disp_tel_parameter"].quantity[rows])
    positions = {s: direction(table[rows], s * norm) for s in (1, -1)}
    best, best_signs = np.inf, None
    for signs in product((1, -1), repeat=len(rows)):
        alt = [positions[s][0][i] for i, s in enumerate(signs)]
        az = [positions[s][1][i] for i, s in enumerate(signs)]
        distance = sum(angular_separation(az[i], alt[i], az[j], alt[j]).to_value(u.deg)
                       for i in range(len(rows)) for j in range(i + 1, len(rows)))
        if distance < best:
            best, best_signs = distance, np.array(signs)
    return best_signs


@pytest.mark.parametrize("n_tel", [2, 3])
def test_combiner_table(n_tel):
    """weighted mean of the closest positions; events with one telescope invalid (min_telescopes 2)"""
    table = mono_table(n_tel=n_tel)
    stereo = StereoDispSignClosestCombiner(prefix="disp", property=GEOMETRY).predict_table(table)

    corrected = table.copy()
    n_valid = {}
    for event_id in np.unique(table["event_id"]):
        rows = np.flatnonzero((table["event_id"] == event_id) & table["disp_tel_is_valid"])
        n_valid[event_id] = len(rows)
        if len(rows) >= 2:
            disp = np.abs(table["disp_tel_parameter"].quantity[rows])
            corrected["disp_tel_parameter"][rows] = (brute_force_signs(table, rows) * disp).to_value(u.deg)
    corrected["disp_tel_alt"], corrected["disp_tel_az"] = direction(corrected, corrected["disp_tel_parameter"].quantity)
    expected = StereoMeanCombiner(prefix="disp", property=GEOMETRY, weights="intensity").predict_table(corrected)

    two = np.array([n_valid[e] >= 2 for e in stereo["event_id"]])
    np.testing.assert_array_equal(stereo["disp_is_valid"], two)
    assert np.any(two) and not np.all(two)
    separation = angular_separation(stereo["disp_az"].quantity, stereo["disp_alt"].quantity,
                                    expected["disp_az"].quantity, expected["disp_alt"].quantity)[two]
    assert np.all(separation < 1e-6 * u.deg)
    assert np.all(np.isnan(stereo["disp_alt"][~two]))
    # at least one sign changed
    assert np.any(corrected["disp_tel_parameter"] != table["disp_tel_parameter"])
    # the input table is not modified
    np.testing.assert_array_equal(mono_table(n_tel=n_tel)["disp_tel_alt"], table["disp_tel_alt"])


def test_combiner_min_telescopes():
    """min_telescopes 1: the events with one telescope keep the direction of the classifier sign"""
    table = mono_table()
    stereo = StereoDispSignClosestCombiner(prefix="disp", property=GEOMETRY, min_telescopes=1).predict_table(table)
    mean = StereoMeanCombiner(prefix="disp", property=GEOMETRY, weights="intensity").predict_table(table)
    np.testing.assert_array_equal(stereo["disp_is_valid"], mean["disp_is_valid"])
    one = np.array([np.count_nonzero((table["event_id"] == e) & table["disp_tel_is_valid"]) == 1
                    for e in stereo["event_id"]])
    assert np.any(one)
    np.testing.assert_allclose(stereo["disp_alt"][one], mean["disp_alt"][one])


def test_combiner_event():
    """the event loop version gives the direction of the table version"""
    _check_event_loop(StereoDispSignClosestCombiner)


def _check_event_loop(combiner_cls):
    table = mono_table(n_events=1)
    table["disp_tel_is_valid"] = True
    combiner = combiner_cls(prefix="disp", property=GEOMETRY, min_telescopes=2)
    expected = combiner.predict_table(table)[0]

    event = ArrayEventContainer()
    for row in table:
        tel_id = row["tel_id"]
        event.dl1.tel[tel_id].parameters = ImageParametersContainer()
        hillas = event.dl1.tel[tel_id].parameters.hillas
        hillas.fov_lon, hillas.fov_lat = row["hillas_fov_lon"] * u.deg, row["hillas_fov_lat"] * u.deg
        hillas.psi, hillas.intensity = row["hillas_psi"] * u.deg, row["hillas_intensity"]
        event.monitoring.tel[tel_id].pointing.altitude = row["telescope_pointing_altitude"] * u.deg
        event.monitoring.tel[tel_id].pointing.azimuth = row["telescope_pointing_azimuth"] * u.deg
        event.dl2.tel[tel_id].disp["disp"] = DispContainer(parameter=row["disp_tel_parameter"] * u.deg)
        event.dl2.tel[tel_id].geometry["disp"] = ReconstructedGeometryContainer(
            alt=row["disp_tel_alt"] * u.deg, az=row["disp_tel_az"] * u.deg, is_valid=True)
    combiner(event)
    stereo = event.dl2.stereo.geometry["disp"]
    assert stereo.is_valid and stereo.telescopes == [1, 2]
    separation = angular_separation(stereo.az, stereo.alt, expected["disp_az"] * u.deg, expected["disp_alt"] * u.deg)
    assert separation < 1e-6 * u.deg
    # the directions of the telescopes keep their predicted signs
    for row in table:
        assert event.dl2.tel[row["tel_id"]].geometry["disp"].alt == row["disp_tel_alt"] * u.deg

    # one valid telescope: invalid stereo direction
    event.dl2.tel[2].geometry["disp"].is_valid = False
    combiner(event)
    assert not event.dl2.stereo.geometry["disp"].is_valid


def test_axes_crossing_disp_signs():
    """each telescope takes the sign of its disp whose position is the closest to the crossing of the axes"""
    # major axes crossing at (0, 0): tel 1 centroid (1, 0) along lon, tel 2 centroid (0, 2) along lat
    hillas = [dict(fov_lon=1 * u.deg, fov_lat=0 * u.deg, psi=0 * u.deg),
              dict(fov_lon=0 * u.deg, fov_lat=2 * u.deg, psi=90 * u.deg)]
    for disps in ([0.9, 1.5], [-0.9, 1.5], [-0.9, -1.5]):
        np.testing.assert_array_equal(axes_crossing_disp_signs(hillas, disps * u.deg), [-1, -1])
    hillas[0]["fov_lon"], hillas[1]["fov_lat"] = -1 * u.deg, -2 * u.deg
    np.testing.assert_array_equal(axes_crossing_disp_signs(hillas, [-0.9, -1.5] * u.deg), [1, 1])
    # parallel axes: no crossing, the predicted signs
    parallel = [dict(fov_lon=0, fov_lat=0, psi=30), dict(fov_lon=0, fov_lat=1, psi=30)]
    np.testing.assert_array_equal(axes_crossing_disp_signs(parallel, [0.5, -0.5]), [1, -1])
    with pytest.raises(ValueError, match="two telescopes"):
        axes_crossing_disp_signs(hillas[:1], [0.9])
    with pytest.raises(NotImplementedError, match="more than two telescopes"):
        axes_crossing_disp_signs([*hillas, hillas[0]], [0.9, 1.5, 0.3])


def test_crossing_combiner():
    """former StereoMARSCombiner: signs of the axes crossing, events with one telescope kept"""
    combiner = StereoDispSignCrossingCombiner(prefix="disp", property=GEOMETRY)
    table = mono_table()
    stereo = combiner.predict_table(table)
    corrected = table.copy()
    for event_id in np.unique(table["event_id"]):
        rows = np.flatnonzero((table["event_id"] == event_id) & table["disp_tel_is_valid"])
        if len(rows) == 2:
            hillas = [{name: table[f"hillas_{name}"].quantity[row] for name in ("fov_lon", "fov_lat", "psi")}
                      for row in rows]
            disp = table["disp_tel_parameter"].quantity[rows]
            corrected["disp_tel_parameter"][rows] = (axes_crossing_disp_signs(hillas, disp) * np.abs(disp)).to_value(u.deg)
    corrected["disp_tel_alt"], corrected["disp_tel_az"] = direction(corrected, corrected["disp_tel_parameter"].quantity)
    expected = StereoMeanCombiner(prefix="disp", property=GEOMETRY, weights="intensity").predict_table(corrected)
    np.testing.assert_array_equal(stereo["disp_is_valid"], expected["disp_is_valid"])
    valid = stereo["disp_is_valid"]
    separation = angular_separation(stereo["disp_az"].quantity, stereo["disp_alt"].quantity,
                                    expected["disp_az"].quantity, expected["disp_alt"].quantity)[valid]
    assert np.all(separation < 1e-6 * u.deg)


def test_crossing_combiner_more_than_two_telescopes():
    table = mono_table(n_events=2, n_tel=3)
    table["disp_tel_is_valid"] = True
    with pytest.raises(NotImplementedError):
        StereoDispSignCrossingCombiner(prefix="disp", property=GEOMETRY).predict_table(table)


@pytest.mark.parametrize("combiner_cls", [StereoDispSignClosestCombiner, StereoDispSignCrossingCombiner])
def test_combiners_event_loop(combiner_cls):
    """the event loop gives the direction of the table version, for both combiners"""
    _check_event_loop(combiner_cls)
