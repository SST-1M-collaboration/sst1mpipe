"""
Positions of the telescopes: subarray of sst1mpipe/resources/instrument/sst1m_array.h5, stored in the
files of the pipeline (DL1, DL2, photon lists) and read with SubarrayDescription.from_hdf
"""
import sys

import astropy.units as u
import matplotlib
import numpy as np
import pandas as pd
import pytest
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from ctapipe.instrument import SubarrayDescription

from conftest import OBSERVATION_RUNS, run_process
from sst1mpipe.io import write_photon_list
from sst1mpipe.resources import MC_CONFIG_FILES, SUBARRAY_FILE, TEST_DATA_DIR
from sst1mpipe.scripts import sst1mpipe_obsplan
from sst1mpipe.utils import camera_to_altaz

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


@pytest.mark.parametrize("subarray_file", [None, "dl1"])
def test_obsplan_subarray_file(tmp_path, monkeypatch, subarray_file):
    """sst1mpipe_obsplan uses the site of sst1m_array.h5, or of another subarray file"""
    matplotlib.use("Agg")
    argv = ["sst1mpipe_obsplan", "--date", "20260120", "--ra", "166.1", "--dec", "38.2", "--out-dir", str(tmp_path)]
    if subarray_file == "dl1":
        dl1 = tmp_path / "events.dl1.h5"
        run_process(OBSERVATION_RUNS[22]["path"], dl1, "--max-events=5")
        argv += ["--subarray-file", str(dl1)]
    monkeypatch.setattr(sys, "argv", argv)
    sst1mpipe_obsplan.main()
    assert (tmp_path / "sst1m_observability_2026-01-20.png").exists()
