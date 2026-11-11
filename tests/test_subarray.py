"""
Positions of the telescopes: subarray of sst1mpipe/resources/instrument/sst1m_array.h5, stored in the
files of the pipeline (DL1, DL2, photon lists) and read with SubarrayDescription.from_hdf
"""

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from ctapipe.instrument import SubarrayDescription

from conftest import OBSERVATION_RUNS, run_process
from sst1mpipe.io import write_photon_list
from sst1mpipe.resources import MC_CONFIG_FILES, SUBARRAY_FILE, TEST_DATA_DIR
from sst1mpipe.utils import camera_to_altaz, get_tel_id, merge_subarrays, read_subarray, telescope_location

# longitude (deg), latitude (deg) and height (m) formerly given by the config files
# (telescope_coords, array_center_coords)
CONFIG_TELESCOPES = {21: (14.782924, 49.911800, 510), 22: (14.782078, 49.913084, 510)}
CONFIG_ARRAY_CENTRE = (14.782501, 49.912442, 511.195)
SIMTEL_FILE = TEST_DATA_DIR / "simtel" / "gamma_200_800E3GeV_20_20deg_ATM52_102168.corsika.gz.NSBmed4.simtel.gz"
SUBARRAY = SubarrayDescription.from_hdf(SUBARRAY_FILE)


def geodetic(lon, lat, height):
    return EarthLocation.from_geodetic(lon * u.deg, lat * u.deg, height * u.m)


def distance(location_1, location_2):
    return np.sqrt(sum((a - b) ** 2 for a, b in zip(location_1.to_geocentric(), location_2.to_geocentric(), strict=True)))


def test_subarray_matches_former_config():
    """the positions of the config files were those of the subarray, with the heights rounded to 510 m"""
    for tel_id, coordinates in CONFIG_TELESCOPES.items():
        assert distance(SUBARRAY.tel_earth_locations[tel_id], geodetic(*coordinates)) < 3 * u.m
    assert distance(SUBARRAY.reference_location, geodetic(*CONFIG_ARRAY_CENTRE)) < 0.1 * u.m


def test_ra_dec_unchanged():
    """the RA/Dec of the events do not change with the positions of the subarray (< 0.1 arcsec)"""
    times = Time("2026-01-20T22:00:00") + np.arange(3) * u.min
    alt, az = [60, 65, 70] * u.deg, [170, 180, 190] * u.deg
    for tel_id, coordinates in CONFIG_TELESCOPES.items():
        new = SkyCoord(alt=alt, az=az, frame=AltAz(location=SUBARRAY.tel_earth_locations[tel_id], obstime=times))
        old = SkyCoord(alt=alt, az=az, frame=AltAz(location=geodetic(*coordinates), obstime=times))
        assert np.all(new.icrs.separation(old.icrs) < 0.1 * u.arcsec)


@pytest.mark.parametrize("times", [None, Time("2026-01-20T22:00:00")])
def test_camera_to_altaz(times):
    """the centre of the camera is the pointing direction"""
    location = SUBARRAY.tel_earth_locations[21]
    altaz = camera_to_altaz(0 * u.m, 0 * u.m, 5.6 * u.m, 70 * u.deg, 180 * u.deg, location=location, times=times)
    assert u.isclose(altaz.alt, 70 * u.deg, atol=1e-6 * u.deg)
    assert u.isclose(altaz.az, 180 * u.deg, atol=1e-6 * u.deg)


def test_subarray_of_dl1_file(tmp_path):
    """the DL1 files of sst1mpipe-process store the subarray of the event source"""
    output = tmp_path / "events.dl1.h5"
    run_process(OBSERVATION_RUNS[22]["path"], output, "--max-events=5")
    subarray = SubarrayDescription.from_hdf(output)
    assert list(subarray.tel_ids) == [22]
    assert distance(subarray.tel_earth_locations[22], SUBARRAY.tel_earth_locations[22]) < 1 * u.mm


def test_subarray_of_simulation(tmp_path):
    """the DL1 files of the simulations store the simulated telescopes 1 and 2"""
    output = tmp_path / "gamma.dl1.h5"
    run_process(SIMTEL_FILE, output, "--max-events=5", config=MC_CONFIG_FILES["low"])
    assert list(SubarrayDescription.from_hdf(output).tel_ids) == [1, 2]


def test_photon_list_subarray(tmp_path):
    """the photon lists store the subarray of their DL2 file"""
    output = tmp_path / "events_pl.h5"
    photons = pd.DataFrame({"event_id": [1, 2], "reco_alt": [70.0, 71.0]})
    write_photon_list(photons, output_file=output, telescope="tel_021", subarray=SUBARRAY.select_subarray([21]))
    assert list(SubarrayDescription.from_hdf(output).tel_ids) == [21]


def test_get_tel_id():
    assert get_tel_id("tel_021") == 21 and get_tel_id("tel_001") == 1
    for telescope in ("stereo", "21", "cs1"):
        with pytest.raises(ValueError, match="not a telescope string"):
            get_tel_id(telescope)


def test_telescope_location():
    """the position of a telescope from its string; in stereo, the first telescope of the subarray"""
    assert telescope_location(SUBARRAY, "tel_022") == SUBARRAY.tel_earth_locations[22]
    assert telescope_location(SUBARRAY, "stereo") == SUBARRAY.tel_earth_locations[SUBARRAY.tel_ids[0]]


def test_read_subarray_of_old_file(tmp_path):
    """a file without subarray (written before the subarray was stored) gets the subarray of sst1mpipe"""
    output = tmp_path / "old_events_pl.h5"
    write_photon_list(pd.DataFrame({"event_id": [1], "reco_alt": [70.0]}), output_file=output, telescope="tel_021")
    subarray = read_subarray(output)
    assert list(subarray.tel_ids) == list(SUBARRAY.tel_ids)


def test_merge_subarrays():
    """the subarrays of the mono DL1 files of a stereo file: both telescopes at their positions"""
    merged = merge_subarrays([SUBARRAY.select_subarray([21]), SUBARRAY.select_subarray([22])])
    assert list(merged.tel_ids) == [21, 22]
    for tel_id in (21, 22):
        assert distance(merged.tel_earth_locations[tel_id], SUBARRAY.tel_earth_locations[tel_id]) < 1 * u.mm
    # a telescope in several subarrays keeps its position in the first one
    moved = SubarrayDescription(name="moved", tel_positions={22: [10.0, 0.0, 0.0] * u.m},
                                tel_descriptions={22: SUBARRAY.tel[22]}, reference_location=SUBARRAY.reference_location)
    assert distance(merge_subarrays([moved, SUBARRAY]).tel_earth_locations[22], SUBARRAY.tel_earth_locations[22]) > 1 * u.m


@pytest.mark.parametrize("function", ["apply_models", "stereo_reconstruction"])
def test_subarray_is_required(function):
    reco = pytest.importorskip("sst1mpipe.reco.reco")  # needs pyirf, compatible with the installed astropy
    with pytest.raises(ValueError, match="subarray is required"):
        getattr(reco, function)(pd.DataFrame())
