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
    # about 1.5% of the rows repeat the current_time of the previous row
    assert np.all(np.diff(drive["time"].unix) >= 0)
    for column in STATE_COLUMNS:
        assert drive[column].dtype == bool
    assert drive["is_tracking"].any()


@pytest.mark.parametrize("tel_id", DRIVE_FILES)
def test_drive_radec_mrk421(drives, tel_id):
    """Pointing of the drive while tracking Mrk 421 (wobbles at 0.7 deg) within 2 deg of the source"""
    drive = drives[tel_id]
    mask = (
        drive["is_tracking"]
        & (drive["time"] > Time("2026-01-21T00:30:00"))
        & (drive["time"] < Time("2026-01-21T04:30:00"))
    )
    radec = drive_radec(drive[mask][::100])
    assert np.all(radec.separation(MRK421) < 2 * u.deg)
