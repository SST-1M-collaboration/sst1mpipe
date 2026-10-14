import astropy.units as u
import numpy as np
import pytest
from ctapipe.instrument import SubarrayDescription
from traitlets.config import Config

from conftest import OBSERVATION_RUNS  # noqa: F401 (no IERS download, as for the other tests)
from sst1mpipe.instrument.bending import BendingModel
from sst1mpipe.io import read_drive_logs
from sst1mpipe.resources import SUBARRAY_FILE, TEST_DATA_DIR

COEFFICIENTS = dict(IA=-2.6, IE=0.2, CA=0.1, NPAE=-0.05, AN=0.03, AW=-0.02, TF=0.1, TX=0.01)


@pytest.fixture(scope="module")
def subarray():
    return SubarrayDescription.from_hdf(SUBARRAY_FILE)


def test_no_bending_by_default(subarray):
    model = BendingModel(subarray)
    alt, az = [30, 60, 80] * u.deg, [-250, -100, 10] * u.deg
    sky_alt, sky_az = model.revert(21, alt, az)
    assert u.allclose(sky_alt, alt)
    assert u.allclose(sky_az, az)


def test_revert_inverse_of_apply(subarray):
    model = BendingModel(subarray, config=Config({"BendingModel": COEFFICIENTS}))
    alt, az = np.meshgrid(np.linspace(20, 80, 13), np.linspace(-270, 90, 25))
    alt, az = alt.ravel() * u.deg, az.ravel() * u.deg
    drive_alt, drive_az = model.apply(21, alt, az)
    assert not u.allclose(drive_alt, alt)
    sky_alt, sky_az = model.revert(21, drive_alt, drive_az)
    assert u.allclose(sky_alt, alt, atol=1e-5 * u.deg)
    assert u.allclose(sky_az, az, atol=1e-5 * u.deg)


def test_offsets_per_telescope(subarray):
    config = Config({"BendingModel": {"IA": [["type", "*", 0.0], ["id", 21, -2.6], ["id", 22, -0.6]]}})
    model = BendingModel(subarray, config=config)
    for tel_id, ia in ((21, -2.6), (22, -0.6)):
        d_alt, d_az = model.offsets(tel_id, 60 * u.deg, 100 * u.deg)
        assert u.isclose(d_alt, 0 * u.deg)
        assert u.isclose(d_az, ia * u.deg)


def test_revert_no_convergence(subarray, caplog):
    # large flexure in 1 / tan(alt): no sky position for a drive position close to the horizon
    model = BendingModel(subarray, config=Config({"BendingModel": {"TX": 5.0}}))
    sky_alt, sky_az = model.revert(21, [45, 1] * u.deg, [0, 0] * u.deg)
    assert np.isfinite(sky_alt[0])
    assert u.isclose(model.apply(21, sky_alt[0], sky_az[0])[0], 45 * u.deg, atol=1e-5 * u.deg)
    assert np.isnan(sky_alt[1])
    assert np.isnan(sky_az[1])
    assert "no convergence" in caplog.text


def test_read_drive_logs_bending(subarray):
    path = TEST_DATA_DIR / "drive" / "DriveSystem1_20260120_001.fits"
    model = BendingModel(subarray, config=Config({"BendingModel": {"IA": -2.6, "IE": 0.2}}))
    drive = read_drive_logs(path, subarray=subarray)
    reverted = read_drive_logs(path, subarray=subarray, bending_model=model)
    assert not drive.meta["bending_reverted"]
    assert reverted.meta["bending_reverted"]
    assert u.allclose(drive["pointing"].alt - reverted["pointing"].alt, 0.2 * u.deg, atol=1e-6 * u.deg)
    d_az = (drive["pointing"].az - reverted["pointing"].az).wrap_at(180 * u.deg)
    assert u.allclose(d_az, -2.6 * u.deg, atol=1e-6 * u.deg)
