import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time

from conftest import OBSERVATION_RUNS  # noqa: F401 (no IERS download, as for the other tests)
from sst1mpipe.io import drive_radec, read_drive_logs
from sst1mpipe.io.drive import STATE_COLUMNS
from sst1mpipe.resources import TEST_DATA_DIR

DRIVE_FILES = {
    21: TEST_DATA_DIR / "drive" / "DriveSystem1_20260120_001.fits",
    22: TEST_DATA_DIR / "drive" / "DriveSystem2_20260120_002.fits",
}
MRK421 = SkyCoord(166.1138 * u.deg, 38.2088 * u.deg, frame="icrs")
# largest separation from Mrk 421 (wobbles at 0.7 deg) of the drive position while tracking it,
# which includes the pointing model of the telescope (azimuth offset of about -2.3 deg for tel 21
# and -0.6 deg for tel 22)
MAX_SEPARATION_MRK421 = {21: 1.87 * u.deg, 22: 1.01 * u.deg}


@pytest.fixture(scope="module")
def drives():
    return {tel_id: read_drive_logs(path) for tel_id, path in DRIVE_FILES.items()}


@pytest.mark.parametrize("tel_id", DRIVE_FILES)
def test_read_drive_logs(drives, tel_id):
    drive = drives[tel_id]
    assert drive.meta["tel_id"] == tel_id
    assert drive.meta["frame"] == "fk5"
    assert drive.meta["equinox"] == "J2000"
    assert isinstance(drive["time"], Time)
    assert isinstance(drive["pointing"].frame, AltAz)
    # about 1.5% of the rows repeat the current_time of the previous row, see
    # https://github.com/SST-1M-collaboration/sst1mpipe/issues/211
    assert np.all(np.diff(drive["time"].unix) >= 0)
    for column in STATE_COLUMNS:
        assert drive[column].dtype == bool
    assert drive["is_tracking"].any()


@pytest.mark.parametrize("tel_id", DRIVE_FILES)
def test_drive_radec_mrk421(drives, tel_id):
    """Pointing of the drive while tracking Mrk 421 within MAX_SEPARATION_MRK421 of the source"""
    drive = drives[tel_id]
    mask = (
        drive["is_tracking"]
        & (drive["time"] > Time("2026-01-21T00:30:00"))
        & (drive["time"] < Time("2026-01-21T04:30:00"))
    )
    radec = drive_radec(drive[mask][::100])
    assert np.all(radec.separation(MRK421) < MAX_SEPARATION_MRK421[tel_id])


def test_unknown_telescope():
    """a drive system of a telescope not in the subarray: clear error"""
    from ctapipe.instrument import SubarrayDescription

    from sst1mpipe.resources import SUBARRAY_FILE

    subarray = SubarrayDescription.from_hdf(SUBARRAY_FILE).select_subarray([22])
    with pytest.raises(ValueError, match="telescope 21 .*DriveSystem1.* is not in the subarray"):
        read_drive_logs(DRIVE_FILES[21], subarray=subarray)


def test_tracking_periods_without_tracking(drives):
    """no sample tracking on source: no tracking period"""
    from sst1mpipe.scripts.sst1mpipe_plot_drive_pointing_tool import tracking_periods

    drive = drives[21].copy()
    drive["is_on_source"] = False
    track, periods = tracking_periods(drive, min_tracking=0)
    assert len(track) == 0 and periods == []


def test_sample_durations():
    """the time to the next sample of the period (at most 5 s), the median step for the last one"""
    from sst1mpipe.scripts.sst1mpipe_plot_drive_pointing_tool import sample_durations

    times = np.array([0.0, 1.0, 2.0, 3.0, 100.0, 101.0, 109.0, 500.0])
    periods = [np.arange(4), np.arange(4, 7), np.array([7])]
    np.testing.assert_allclose(sample_durations(times, periods), [1, 1, 1, 1, 1, 5, 3, 0])


def test_field_targets_started_before_the_field():
    """a run started during the slew (before the field) is drawn on it, not one started long before"""
    from sst1mpipe.scripts.sst1mpipe_plot_drive_pointing_tool import field_targets

    start = Time("2026-01-20T20:00:00")
    field = dict(start=start, stop=start + 1 * u.h, center=MRK421)
    targets = [
        dict(name="slew", time=start - 2 * u.min, coord=MRK421),
        dict(name="earlier", time=start - 1 * u.h, coord=MRK421),
        dict(name="other field", time=start + 10 * u.min, coord=MRK421.directional_offset_by(0 * u.deg, 20 * u.deg)),
    ]
    assert [target["name"] for target in field_targets(field, targets, 10 * u.deg)] == ["slew"]


@pytest.mark.parametrize("frame", ["icrs", "fk5"])
def test_wcs_frame(frame):
    """the maps accept sky frames with and without equinox"""
    from sst1mpipe.scripts.sst1mpipe_plot_drive_pointing_tool import field_wcs

    wcs, shape = field_wcs(MRK421.transform_to(frame), half_width=1.0)
    assert wcs.wcs.radesys == frame.upper() and shape == (100, 100)
