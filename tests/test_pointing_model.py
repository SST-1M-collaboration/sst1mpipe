import json

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time
from ctapipe.core import run_tool
from ctapipe.instrument import SubarrayDescription
from traitlets.config import Config

from conftest import OBSERVATION_RUNS  # noqa: F401 (no IERS download, as for the other tests)
from sst1mpipe.instrument.pointing_model import PointingModel, TPointPointingModel
from sst1mpipe.io import read_drive_logs
from sst1mpipe.resources import SST1M_POINTING_MODEL_CONFIG_FILE, SUBARRAY_FILE, TEST_DATA_DIR

# arcsec, of the size of the corrections of the SST-1M drives
COEFFICIENTS = dict(ia=-8360.0, ie=970.0, an=36.0, aw=220.0, ca=-250.0, npae=-180.0, tf=360.0, tx=40.0)
DRIVE_FILES = [TEST_DATA_DIR / "drive" / "DriveSystem1_20260120_001.fits", TEST_DATA_DIR / "drive" / "DriveSystem2_20260120_002.fits"]


@pytest.fixture(scope="module")
def subarray():
    return SubarrayDescription.from_hdf(SUBARRAY_FILE)


def altaz(subarray, alt, az, tel_id=21):
    frame = AltAz(obstime=Time("2026-01-21T03:00:00"), location=subarray.tel_earth_locations[tel_id])
    return SkyCoord(alt=u.Quantity(alt, u.deg), az=u.Quantity(az, u.deg), frame=frame)


def tpoint(subarray, **coefficients):
    return TPointPointingModel(subarray, config=Config({"TPointPointingModel": coefficients}))


def test_pointing_model_is_abstract(subarray):
    with pytest.raises(TypeError):
        PointingModel(subarray)
    assert isinstance(PointingModel.from_name("TPointPointingModel", subarray=subarray), TPointPointingModel)


def test_no_correction_by_default(subarray):
    model = TPointPointingModel(subarray)
    pointing = altaz(subarray, [30, 60, 80], [-250, -100, 10])
    for coord in (model.revert(21, pointing), model.apply(21, pointing)):
        assert u.allclose(coord.alt, pointing.alt)
        assert u.allclose(coord.az, pointing.az)


@pytest.mark.parametrize(
    "term, alt, az, d_alt, d_az",
    [
        # TPOINT terms (pytelpoint.transform.azel_model)
        ("ia", 45, 30, 0, -1),
        ("ie", 45, 30, 1, 0),
        ("an", 45, 0, -1, 0),
        ("an", 45, 90, 0, -1),
        ("aw", 45, 90, 1, 0),
        ("aw", 45, 0, 0, -1),
        ("ca", 60, 30, 0, -2),
        ("npae", 45, 30, 0, -1),
        ("tf", 60, 30, -0.5, 0),
        ("tx", 45, 30, -1, 0),
        # centring terms (signs of the SST-1M drives)
        ("aces", 45, 90, 0, -1),
        ("acec", 45, 0, 0, 1),
        ("eces", 30, 30, 0.5, 0),
        ("ecec", 60, 30, 0.5, 0),
    ],
)
def test_tpoint_terms(subarray, term, alt, az, d_alt, d_az):
    model = tpoint(subarray, **{term: 1.0})
    corrections = model.corrections(21, alt * u.deg, az * u.deg)
    assert u.isclose(corrections[0], d_alt * u.arcsec, atol=1e-9 * u.arcsec)
    assert u.isclose(corrections[1], d_az * u.arcsec, atol=1e-9 * u.arcsec)
    # sky = drive + corrections(drive)
    sky = model.revert(21, altaz(subarray, alt, az))
    assert u.isclose(sky.alt, alt * u.deg + d_alt * u.arcsec)


def test_apply_inverse_of_revert(subarray, caplog):
    model = tpoint(subarray, **COEFFICIENTS)
    alt, az = np.meshgrid(np.linspace(20, 75, 12), np.linspace(0, 355, 72))
    drive = altaz(subarray, alt.ravel(), az.ravel())
    sky = model.revert(21, drive)
    assert np.all(sky.separation(drive) > 0.3 * u.deg)
    # first order inverse: within the tolerance up to 75 deg of altitude for these corrections
    assert np.all(model.apply(21, sky).separation(drive) < model.tolerance * u.arcsec)
    assert "error above" not in caplog.text


def test_apply_tolerance(subarray, caplog):
    # large tilt of the azimuth axis close to the zenith: the first order inverse is not precise enough
    model = tpoint(subarray, an=0.3 * 3600)
    sky = altaz(subarray, [45, 85], [30, 30])
    drive = model.apply(21, sky)
    assert np.all(np.isfinite(drive.alt))
    distance = model.revert(21, drive).separation(sky)
    assert distance[0] < model.tolerance * u.arcsec
    assert distance[1] > model.tolerance * u.arcsec
    assert "Pointing model of telescope 21: 1 positions inverted with an error above 30.0 arcsec" in caplog.text


def test_coefficients_per_telescope(subarray):
    model = tpoint(subarray, ia=[["type", "*", 0.0], ["id", 21, 3600.0], ["id", 22, 1800.0]])
    for tel_id, d_az in ((21, -1.0), (22, -0.5)):
        drive = altaz(subarray, 60, 100, tel_id)
        sky = model.revert(tel_id, drive)
        assert u.isclose(sky.alt, drive.alt)
        assert u.isclose((sky.az - drive.az).wrap_at(180 * u.deg), d_az * u.deg)


def test_revert_drive_logs(subarray):
    model = tpoint(subarray, ia=-3600.0, ie=720.0)
    drive = read_drive_logs(DRIVE_FILES[0], subarray=subarray)
    sky = model.revert(21, drive["pointing"])
    assert sky.frame.is_equivalent_frame(drive["pointing"].frame)
    tracking = drive["is_tracking"]
    assert u.allclose(sky.alt[tracking] - drive["pointing"].alt[tracking], 0.2 * u.deg, atol=1e-6 * u.deg)
    d_az = (sky.az - drive["pointing"].az).wrap_at(180 * u.deg)
    assert u.allclose(d_az[tracking], 1.0 * u.deg, atol=1e-6 * u.deg)


@pytest.mark.parametrize("revert_pointing_model", [False, True])
def test_plot_drive_pointing_tool(tmp_path, revert_pointing_model):
    from sst1mpipe.scripts.sst1mpipe_plot_drive_pointing_tool import PlotDrivePointingTool

    argv = ["-i", str(DRIVE_FILES[0]), "-i", str(DRIVE_FILES[1]), "-o", str(tmp_path / "plots")]
    if revert_pointing_model:
        config = tmp_path / "pointing_model.json"
        config.write_text(json.dumps({"TPointPointingModel": {"ia": [["type", "*", 0], ["id", 21, -8360], ["id", 22, -2070]]}}))
        argv += ["--revert-pointing-model", "--config", str(config)]
    assert run_tool(PlotDrivePointingTool(), argv=argv, cwd=tmp_path) == 0
    for name in ("drive_radec_vs_time.png", "drive_radec_map.png"):
        assert (tmp_path / "plots" / name).stat().st_size > 0


@pytest.mark.parametrize("tel_id, max_separation", [(21, 1 * u.arcmin), (22, 2 * u.arcmin)])
def test_sst1m_pointing_model(subarray, tel_id, max_separation):
    """
    Coefficients of the SST-1M drives (SST1M_POINTING_MODEL_CONFIG_FILE): pointing on the Mrk 421
    wobble W1 (TARGET of the zfits runs) within max_separation
    """
    model = TPointPointingModel(subarray, config=Config(json.loads(SST1M_POINTING_MODEL_CONFIG_FILE.read_text())))
    assert model.ia.tel[tel_id] != 0
    drive = read_drive_logs(DRIVE_FILES[tel_id - 21], subarray=subarray)
    tracking = (
        drive["is_tracking"] & drive["is_on_source"]
        & (drive["time"] > Time("2026-01-21T00:30:00")) & (drive["time"] < Time("2026-01-21T04:30:00"))
    )
    sky = model.revert(tel_id, drive["pointing"][tracking]).transform_to("icrs")
    w1 = SkyCoord(166.9948 * u.deg, 38.1053 * u.deg, frame="icrs")
    on_w1 = sky.ra > 166.5 * u.deg
    assert np.count_nonzero(on_w1) > 1000
    assert np.all(SkyCoord(sky.ra, sky.dec)[on_w1].separation(w1) < max_separation)
