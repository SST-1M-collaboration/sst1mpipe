import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time
from ctapipe.core import run_tool
from ctapipe.instrument import SubarrayDescription
from traitlets.config import Config

from conftest import OBSERVATION_RUNS  # noqa: F401 (no IERS download, as for the other tests)
from sst1mpipe.instrument.bending import BendingModel, MeeksBendingModel
from sst1mpipe.io import read_drive_logs
from sst1mpipe.resources import SUBARRAY_FILE, TEST_DATA_DIR

# deg, of the size of the pointing errors of the SST-1M drives
PARAMETERS = dict(
    azimuth_axis_tilt=0.05, level_axis_azimuth=40.0, elevation_axis_skew=-0.05, collimation=0.07,
    elevation_coefficient=1e-3, azimuth_bias=-2.3, elevation_bias=0.27,
)
DRIVE_FILES = [TEST_DATA_DIR / "drive" / "DriveSystem1_20260120_001.fits", TEST_DATA_DIR / "drive" / "DriveSystem2_20260120_002.fits"]


@pytest.fixture(scope="module")
def subarray():
    return SubarrayDescription.from_hdf(SUBARRAY_FILE)


def altaz(subarray, alt, az, tel_id=21):
    frame = AltAz(obstime=Time("2026-01-21T03:00:00"), location=subarray.tel_earth_locations[tel_id])
    return SkyCoord(alt=u.Quantity(alt, u.deg), az=u.Quantity(az, u.deg), frame=frame)


def meeks(subarray, **parameters):
    return MeeksBendingModel(subarray, config=Config({"MeeksBendingModel": parameters}))


def test_bending_model_is_abstract(subarray):
    with pytest.raises(TypeError):
        BendingModel(subarray)
    assert isinstance(BendingModel.from_name("MeeksBendingModel", subarray=subarray), MeeksBendingModel)


def test_no_pointing_error_by_default(subarray):
    model = MeeksBendingModel(subarray)
    pointing = altaz(subarray, [30, 60, 80], [-250, -100, 10])
    for coord in (model.apply(21, pointing), model.revert(21, pointing)):
        assert u.allclose(coord.alt, pointing.alt)
        assert u.allclose(coord.az, pointing.az)


@pytest.mark.parametrize(
    "parameters, alt, az, d_alt, d_az",
    [
        # azimuth-axis tilt, eq. 3: tilted toward kappa_A + 90 deg
        (dict(azimuth_axis_tilt=1.0), 45, 0, 0, -1),
        (dict(azimuth_axis_tilt=1.0), 45, 90, 1, 0),
        (dict(azimuth_axis_tilt=1.0, level_axis_azimuth=90.0), 45, 90, 0, -1),
        (dict(azimuth_axis_tilt=1.0, level_axis_azimuth=90.0), 45, 180, 1, 0),
        # elevation-axis skew, eq. 5
        (dict(elevation_axis_skew=1.0), 45, 30, 0, 1),
        # collimation error, eq. 6
        (dict(collimation=1.0), 60, 30, 0, -2),
        # gravitational deflection, eq. 7, h in deg
        (dict(elevation_coefficient=1e-3), 60, 30, 0.06, 0),
        # biases of the encoders, eqs. 8 and 9
        (dict(azimuth_bias=1.0), 45, 30, 0, 1),
        (dict(elevation_bias=1.0), 45, 30, 1, 0),
    ],
)
def test_meeks_pointing_errors(subarray, parameters, alt, az, d_alt, d_az):
    d = meeks(subarray, **parameters).pointing_errors(21, alt * u.deg, az * u.deg)
    assert u.isclose(d[0], d_alt * u.deg, atol=1e-12 * u.deg)
    assert u.isclose(d[1], d_az * u.deg, atol=1e-12 * u.deg)


def test_revert_inverse_of_apply(subarray):
    model = meeks(subarray, **PARAMETERS)
    alt, az = np.meshgrid(np.linspace(20, 80, 13), np.linspace(0, 355, 72))
    sky = altaz(subarray, alt.ravel(), az.ravel())
    drive = model.apply(21, sky)
    assert np.all(drive.separation(sky) > 0.3 * u.deg)
    assert np.all(model.revert(21, drive).separation(sky) < 1e-2 * u.arcsec)


def test_parameters_per_telescope(subarray):
    model = meeks(subarray, azimuth_bias=[["type", "*", 0.0], ["id", 21, 1.0], ["id", 22, 0.5]])
    for tel_id, d_az in ((21, 1.0), (22, 0.5)):
        sky = altaz(subarray, 60, 100, tel_id)
        drive = model.apply(tel_id, sky)
        assert u.isclose(drive.alt, sky.alt)
        assert u.isclose((drive.az - sky.az).wrap_at(180 * u.deg), d_az * u.deg)


def test_revert_no_convergence(subarray, caplog):
    model = meeks(subarray, collimation=5.0, azimuth_axis_tilt=2.0)
    sky = model.revert(21, altaz(subarray, [45, 89.9], [30, 30]))
    assert np.isfinite(sky.alt[0])
    assert model.apply(21, sky[0]).separation(altaz(subarray, 45, 30)) < 1e-2 * u.arcsec
    assert np.isnan(sky.alt[1])
    assert "no convergence" in caplog.text


def test_revert_drive_logs(subarray):
    model = meeks(subarray, azimuth_bias=1.0, elevation_bias=0.2)
    drive = read_drive_logs(DRIVE_FILES[0], subarray=subarray)
    sky = model.revert(21, drive["pointing"])
    assert sky.frame.is_equivalent_frame(drive["pointing"].frame)
    tracking = drive["is_tracking"]
    assert u.allclose(sky.alt[tracking] - drive["pointing"].alt[tracking], -0.2 * u.deg, atol=1e-6 * u.deg)
    d_az = (sky.az - drive["pointing"].az).wrap_at(180 * u.deg)
    assert u.allclose(d_az[tracking], -1.0 * u.deg, atol=1e-6 * u.deg)


@pytest.mark.parametrize("revert_bending", [False, True])
def test_plot_drive_pointing_tool(tmp_path, revert_bending):
    from sst1mpipe.scripts.sst1mpipe_plot_drive_pointing_tool import PlotDrivePointingTool

    argv = ["-i", str(DRIVE_FILES[0]), "-i", str(DRIVE_FILES[1]), "-o", str(tmp_path / "plots")]
    if revert_bending:
        config = tmp_path / "bending.json"
        config.write_text('{"MeeksBendingModel": {"azimuth_bias": [["type", "*", 0], ["id", 21, -2.3], ["id", 22, -0.6]]}}')
        argv += ["--revert-bending", "--config", str(config)]
    assert run_tool(PlotDrivePointingTool(), argv=argv, cwd=tmp_path) == 0
    for name in ("drive_radec_vs_time.png", "drive_radec_map.png"):
        assert (tmp_path / "plots" / name).stat().st_size > 0
