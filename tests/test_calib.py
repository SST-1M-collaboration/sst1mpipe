
import numpy as np
import pytest
import json
from copy import deepcopy

from ctapipe.calib import CameraCalibrator
from astropy.table import Table
from astropy.time import Time
from ctapipe.containers import (
    CameraCalibrationContainer,
    PixelStatus,
    ArrayEventContainer,
    R0CameraContainer,
    R1CameraContainer,
    SimulatedCameraContainer,
    SimulatedEventContainer,
)
from ctapipe.instrument import SubarrayDescription
from traitlets.config import Config

from sst1mpipe.calib import ImageSaturationCorrector, R0R1Calibrator
from sst1mpipe.calib.calib import DEFAULT_CALIBRATION_FILES, DEFAULT_WINDOW_FILES
from conftest import DARK_FILE
from sst1mpipe.io import load_config, translate_legacy_calibration_config
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.monitoring import CameraCoefficientsCalculator
from sst1mpipe.utils import get_subarray
from sst1mpipe.resources import (
    CALIBRATION_DIR,
    CONFIG_DIR,
    DATA_CONFIG_FILE,
    MC_CONFIG_FILES,
    MONITORING_CONFIG_FILE,
    PDE_CORRECTION_FACTORS_FILE,
    WINDOW_DIR,
)

CONFIG = load_config(DATA_CONFIG_FILE, ismc=False)
TEL_ID = 22
CALIBRATION_FILE_TEL_2 = str((CALIBRATION_DIR / DEFAULT_CALIBRATION_FILES[22]))


@pytest.fixture(scope="module")
def event():
    source = SST1MEventSource(input_url=DARK_FILE, max_events=1)
    return next(iter(source))


def calibrator(**settings):
    return R0R1Calibrator(subarray=get_subarray(), config=Config({"R0R1Calibrator": settings}))


def without_coefficients(event):
    """the event without camera calibration coefficients (no monitoring source)"""
    event.monitoring.tel[TEL_ID].camera.coefficients = CameraCalibrationContainer()
    return event


def coefficients_calculator(**settings):
    return CameraCoefficientsCalculator(subarray=get_subarray(), config=Config({"CameraCoefficientsCalculator": settings}))


def waveform_statistics(pedestal_std, outlier_mask=None):
    """statistics of the ADC samples of one chunk of pedestal events (as compute_pixel_statistics)"""
    std = np.asarray(pedestal_std, dtype=np.float64).reshape(1, 1, -1)
    return Table(dict(
        time_start=Time(["2026-01-20T22:00:00"]),
        time_end=Time(["2026-01-20T22:00:10"]),
        std=std,
        outlier_mask=np.zeros(std.shape, dtype=bool) if outlier_mask is None else np.reshape(outlier_mask, std.shape),
        is_valid=[True],
    ))


def test_r0_r1_dl1_calibration():

    source = SST1MEventSource(input_url=DARK_FILE, max_events=3)
    n_pixels = source.subarray.tel[TEL_ID].camera.readout.n_pixels
    calibrator_r0_r1 = R0R1Calibrator(subarray=source.subarray, config=CONFIG)
    r1_dl1_calibrator = CameraCalibrator(subarray=source.subarray, config=CONFIG)

    for event in source:
        r0_waveform = event.r0.tel[TEL_ID].waveform.copy()

        calibrator_r0_r1(event, TEL_ID)
        r1 = event.r1.tel[TEL_ID]

        # ctapipe shape: (n_channels, n_pixels, n_samples)
        assert r1.waveform.shape == r0_waveform.shape == (1, n_pixels, r0_waveform.shape[-1])
        bad_pixels = event.monitoring.tel[TEL_ID].camera.coefficients.outlier_mask[0]
        assert bad_pixels.any()
        assert calibrator_r0_r1.n_bad_pixels[TEL_ID] == bad_pixels.sum()
        # the raw data are not modified by the calibration
        np.testing.assert_array_equal(event.r0.tel[TEL_ID].waveform, r0_waveform)

        r1_dl1_calibrator(event)
        assert event.dl1.tel[TEL_ID].image.shape == (n_pixels, )
        # no saturated pixels in the dark run of the test file
        assert not ImageSaturationCorrector(subarray=source.subarray, config=CONFIG)(event, TEL_ID)


def test_all_telescopes_with_r0_data_are_calibrated(event):

    event.r1.tel.clear()
    calibrator()(event)

    assert list(event.r1.tel.keys()) == [TEL_ID]


def test_pedestal_subtraction_and_dc_to_pe(event):

    calibrator_r0_r1 = calibrator(flag_bad_calibration_pixels=False)
    calibrator_r0_r1(without_coefficients(event), TEL_ID)

    r0 = event.r0.tel[TEL_ID]
    dc_to_pe, _ = calibrator_r0_r1.calibration_parameters(TEL_ID)
    window = calibrator_r0_r1.window_transmittance(TEL_ID)
    expected = (r0.waveform[0] - r0.pedestal[:, np.newaxis]) / dc_to_pe[:, np.newaxis] / window[:, np.newaxis]
    np.testing.assert_allclose(event.r1.tel[TEL_ID].waveform[0], expected)
    assert calibrator_r0_r1.n_bad_pixels[TEL_ID] == 0


def test_no_voltage_drop_correction_without_coefficients(event):
    """without coefficients of the monitoring, the factor of the coefficients is not filled"""
    calibrator()(without_coefficients(event), TEL_ID)

    assert event.monitoring.tel[TEL_ID].camera.coefficients.factor is None


@pytest.mark.parametrize("flag_bad_calibration", [False, True])
def test_bad_pixels(event, flag_bad_calibration):
    """the pixels with bad calibration parameters are added to the outliers of the monitoring (dead pixels)"""
    dead_pixel = 100
    dead = np.zeros((1, 1296), dtype=bool)
    dead[0, dead_pixel] = True
    event = without_coefficients(event)
    event.monitoring.tel[TEL_ID].camera.coefficients = CameraCalibrationContainer(outlier_mask=dead)
    calibrator_r0_r1 = calibrator(flag_bad_calibration_pixels=flag_bad_calibration)
    calibrator_r0_r1(event, TEL_ID)

    _, mask_bad_calibration = calibrator_r0_r1.calibration_parameters(TEL_ID)
    flagged = event.monitoring.tel[TEL_ID].camera.coefficients.outlier_mask[0]
    assert mask_bad_calibration.any() and not mask_bad_calibration[dead_pixel]
    assert np.all(flagged[mask_bad_calibration] == flag_bad_calibration)
    assert flagged[dead_pixel]
    assert flagged.sum() == 1 + flag_bad_calibration * mask_bad_calibration.sum()
    # the outliers of the monitoring are not modified
    assert dead.sum() == 1


@pytest.mark.parametrize("correction", ["none", "global", "pixelwise"])
def test_voltage_drop_coefficients(correction):

    pedestal_std = np.linspace(3, 6, 1296)
    coefficients = coefficients_calculator(voltage_drop_correction=correction)(waveform_statistics(pedestal_std), TEL_ID)
    # the charges are multiplied by the factor: 1 / charge drop
    factor = np.asarray(coefficients["factor"][0, 0])
    if correction == "none":
        assert np.all(factor == 1)
    elif correction == "global":
        # the charges are increased, more for a higher NSB (pedestal variance)
        assert np.all(factor == factor[0]) and factor[0] > 1
    else:
        assert np.all(np.diff(factor) > 0)
    assert not np.any(coefficients["outlier_mask"])
    assert coefficients["time"][0] == coefficients["time_end"][0]
    np.testing.assert_array_equal(coefficients["pedestal_offset"], 0)
    np.testing.assert_array_equal(coefficients["time_shift"], 0)


@pytest.mark.parametrize("flag_dead", [False, True])
def test_dead_pixels(flag_dead):
    """the pixels without statistics (outliers of the monitoring) are dead pixels, a pixel with a low std is not"""
    dead_pixel, low_std_pixel = 100, 101
    pedestal_std = np.full(1296, 4.0)
    pedestal_std[low_std_pixel] = 1.0
    outliers = np.zeros(1296, dtype=bool)
    outliers[dead_pixel] = True
    statistics = waveform_statistics(pedestal_std, outliers)
    coefficients = coefficients_calculator(flag_dead_pixels=flag_dead)(statistics, TEL_ID)

    flagged = coefficients["outlier_mask"][0, 0]
    assert flagged[dead_pixel] == flag_dead
    assert flagged.sum() == flag_dead


def test_pixels_without_statistics():
    """the pixels without statistics are not used for the global voltage drop correction"""
    outliers = np.zeros(1296, dtype=bool)
    outliers[:10] = True
    pedestal_std = np.full(1296, 4.0)
    pedestal_std[:10] = 100
    calculator = coefficients_calculator(voltage_drop_correction="global")

    factor = calculator(waveform_statistics(pedestal_std, outliers), TEL_ID)["factor"]
    np.testing.assert_allclose(factor, calculator(waveform_statistics(np.full(1296, 4.0)), TEL_ID)["factor"])
    assert np.all(np.isfinite(factor))


def test_calibration_file_per_telescope():

    default = calibrator()
    assert default.calibration_file_path(21).name == DEFAULT_CALIBRATION_FILES[21]
    assert default.calibration_file_path(22).name == DEFAULT_CALIBRATION_FILES[22]

    # use the calibration file of tel 2 for tel 1
    custom = calibrator(calibration_file=[["type", "*", None], ["id", 21, CALIBRATION_FILE_TEL_2]])
    assert str(custom.calibration_file_path(21)) == CALIBRATION_FILE_TEL_2
    assert custom.calibration_file_path(22).name == DEFAULT_CALIBRATION_FILES[22]
    np.testing.assert_array_equal(
        custom.calibration_parameters(21)[0], default.calibration_parameters(22)[0],
    )


def test_wrong_settings():

    with pytest.raises(Exception, match="(?i)trait|exist"):
        calibrator(calibration_file="/does/not/exist.h5")
    with pytest.raises(Exception, match="(?i)trait"):
        coefficients_calculator(voltage_drop_correction="wrong")


@pytest.mark.parametrize("pixelwise, global_, expected", [
    (False, True, "global"),
    (True, False, "pixelwise"),
    (True, True, "pixelwise"),
    (False, False, "none"),
])
def test_translate_legacy_calibration_config(pixelwise, global_, expected):

    legacy = {
        "telescope_calibration": {
            "tel_021": None,
            "tel_022": CALIBRATION_FILE_TEL_2,
            "bad_calib_px_interpolation": True,
            "dynamic_dead_px_interpolation": False,
        },
        "NsbCalibrator": {
            "apply_pixelwise_Vdrop_correction": pixelwise,
            "apply_global_Vdrop_correction": global_,
            "intensity_correction": {"tel_021": 1.0},
            "mc_correction_for_PDE": False,
        },
    }
    config = translate_legacy_calibration_config(legacy)

    assert "telescope_calibration" not in config
    assert config["NsbCalibrator"] == {"intensity_correction": {"tel_021": 1.0}}
    assert config["R0R1Calibrator"] == {
        "calibration_file": [["type", "*", None], ["id", 22, CALIBRATION_FILE_TEL_2]],
        "flag_bad_calibration_pixels": True,
    }
    assert config["CameraCoefficientsCalculator"] == {"flag_dead_pixels": False, "voltage_drop_correction": expected}
    calibrator_r0_r1 = R0R1Calibrator(subarray=get_subarray(), config=Config(config))
    assert str(calibrator_r0_r1.calibration_file_path(22)) == CALIBRATION_FILE_TEL_2
    calculator = CameraCoefficientsCalculator(subarray=get_subarray(), config=Config(config))
    assert calculator.voltage_drop_correction.tel[21] == expected


def test_translate_legacy_mc_pde_correction():

    with pytest.raises(ValueError, match="pde_drop_factor"):
        translate_legacy_calibration_config({"NsbCalibrator": {"mc_correction_for_PDE": True}})


def test_config_without_legacy_settings_is_unchanged():

    config = {"R0R1Calibrator": {"flag_bad_calibration_pixels": False}, "NsbCalibrator": {}}

    assert translate_legacy_calibration_config(config) == config


def test_default_config_settings():

    calibrator_r0_r1 = R0R1Calibrator(subarray=get_subarray(), config=CONFIG)
    calculator = CameraCoefficientsCalculator(subarray=get_subarray(), config=load_config(MONITORING_CONFIG_FILE))

    for tel_id in (21, 22):
        assert calibrator_r0_r1.flag_bad_calibration_pixels.tel[tel_id]
        assert calibrator_r0_r1.calibration_file_path(tel_id).name == DEFAULT_CALIBRATION_FILES[tel_id]
        assert calculator.voltage_drop_correction.tel[tel_id] == "global"
        assert calculator.flag_dead_pixels.tel[tel_id]


# PDE drop correction of the simulations
# PDE files of the low and high NSB simulations, see mc_pde_correction_factors.json
PDE_FILES = {
    "low": {1: "qe_SST1M_5477_ave_TEL1_NSB136.0", 2: "qe_SST1M_5477_ave_TEL2_NSB177.0"},
    "high": {1: "qe_SST1M_5477_ave_TEL1_NSB251.0", 2: "qe_SST1M_5477_ave_TEL2_NSB300.0"},
}


@pytest.fixture(scope="module")
def mc_subarray():
    """SST-1M subarray with the telescope ids of the simulations (1, 2)"""
    subarray = get_subarray()
    return SubarrayDescription(
        "SST1M_MC",
        tel_positions={1: subarray.positions[21], 2: subarray.positions[22]},
        tel_descriptions={1: subarray.tel[21], 2: subarray.tel[22]},
        reference_location=subarray.reference_location,
    )


def simulated_event(tel_ids=(1, 2)):
    event = ArrayEventContainer()
    event.simulation = SimulatedEventContainer()
    for tel_id in tel_ids:
        event.r1.tel[tel_id] = R1CameraContainer(waveform=np.ones((1, 1296, 50), dtype=np.float32))
    return event


def mc_calibrator(mc_subarray, **settings):
    return R0R1Calibrator(subarray=mc_subarray, config=Config({"R0R1Calibrator": settings}))


def test_pde_drop_correction(mc_subarray):

    event = simulated_event()
    calibrator_r0_r1 = mc_calibrator(mc_subarray, pde_drop_factor=[["id", 1, 0.5], ["id", 2, 0.8]])
    calibrator_r0_r1(event)

    # the charges are divided by the PDE drop: factor of the coefficients, the R1 waveforms are kept
    np.testing.assert_allclose(event.monitoring.tel[1].camera.coefficients.factor, 2)
    np.testing.assert_allclose(event.monitoring.tel[2].camera.coefficients.factor, 1.25)
    np.testing.assert_array_equal(event.r1.tel[1].waveform, 1)
    # the R0 data of the simulations are not used
    assert len(event.r0.tel) == 0


def test_pde_drop_correction_single_telescope(mc_subarray):

    event = simulated_event()
    mc_calibrator(mc_subarray, pde_drop_factor=0.5)(event, 2)

    assert event.monitoring.tel[1].camera.coefficients.factor is None
    np.testing.assert_allclose(event.monitoring.tel[2].camera.coefficients.factor, 2)


@pytest.mark.parametrize("settings", [{}, {"pde_drop_factor": None}, {"pde_drop_factor": [["id", 1, 0.5]]}])
def test_no_pde_drop_correction(mc_subarray, settings):

    event = simulated_event(tel_ids=[2])
    calibrator_r0_r1 = mc_calibrator(mc_subarray, **settings)
    calibrator_r0_r1(event)

    assert calibrator_r0_r1.pde_drop(2) is None
    assert event.monitoring.tel[2].camera.coefficients.factor is None
    np.testing.assert_array_equal(event.r1.tel[2].waveform, 1)


@pytest.mark.parametrize("nsb", ["low", "high"])
def test_mc_configs_pde_drop_factors(mc_subarray, nsb):

    with open(PDE_CORRECTION_FACTORS_FILE) as f:
        factors = json.load(f)["mc_correction_for_PDE"]
    config = load_config(MC_CONFIG_FILES[nsb], ismc=True)
    calibrator_r0_r1 = R0R1Calibrator(subarray=mc_subarray, config=config)

    for tel_id in (1, 2):
        expected = factors[f"tel_00{tel_id}"][PDE_FILES[nsb][tel_id]]
        assert calibrator_r0_r1.pde_drop(tel_id) == expected
    assert "intensity_correction" in config["NsbCalibrator"]  # used by sst1mpipe_dl1_dl2


def test_default_mc_config_is_low_nsb():

    # one MC config per NSB level, the low NSB one is used when no config is given
    assert load_config(None, ismc=True) == load_config(MC_CONFIG_FILES["low"])
    assert load_config(None, ismc=True) != load_config(MC_CONFIG_FILES["high"])
    assert not (CONFIG_DIR / 'sst1mpipe_mc_config.json').is_file()


@pytest.mark.parametrize("config_file", [DATA_CONFIG_FILE, *MC_CONFIG_FILES.values()])
def test_no_pde_drop_correction_of_real_telescopes(config_file):

    calibrator_r0_r1 = R0R1Calibrator(subarray=get_subarray(), config=load_config(config_file))

    assert calibrator_r0_r1.pde_drop(21) is None
    assert calibrator_r0_r1.pde_drop(22) is None


def test_pde_drop_correction_of_data(event):

    calibrator(pde_drop_factor=0.5)(without_coefficients(event), TEL_ID)

    np.testing.assert_allclose(event.monitoring.tel[TEL_ID].camera.coefficients.factor, 2)


# bad pixels of the config
BAD_PIXELS = [3, 100, 1201]


def test_static_bad_pixels_of_data(event):

    calibrator_r0_r1 = calibrator(
        flag_bad_calibration_pixels=False, bad_pixels=[["type", "*", []], ["id", TEL_ID, BAD_PIXELS]],
    )
    calibrator_r0_r1(without_coefficients(event), TEL_ID)

    flagged = event.monitoring.tel[TEL_ID].camera.coefficients.outlier_mask[0]
    assert np.flatnonzero(flagged).tolist() == BAD_PIXELS
    assert calibrator_r0_r1.n_bad_pixels[TEL_ID] == len(BAD_PIXELS)


def test_static_bad_pixels_added_to_the_other_bad_pixels(event):

    calibrator()(without_coefficients(event), TEL_ID)
    bad_calibration = event.monitoring.tel[TEL_ID].camera.coefficients.outlier_mask[0].copy()

    calibrator(bad_pixels=[["id", TEL_ID, BAD_PIXELS]])(without_coefficients(event), TEL_ID)
    flagged = event.monitoring.tel[TEL_ID].camera.coefficients.outlier_mask[0]

    assert bad_calibration.any()
    np.testing.assert_array_equal(np.flatnonzero(flagged), np.union1d(np.flatnonzero(bad_calibration), BAD_PIXELS))


def test_static_bad_pixels_of_simulations(mc_subarray):

    event = simulated_event()
    for tel_id in (1, 2):
        event.simulation.tel[tel_id] = SimulatedCameraContainer(true_image=np.ones(1296, dtype=np.int32))
    calibrator_r0_r1 = mc_calibrator(mc_subarray, bad_pixels=[["type", "*", []], ["id", 1, BAD_PIXELS]])
    calibrator_r0_r1(event)

    flagged = event.monitoring.tel[1].camera.coefficients.outlier_mask[0]
    assert np.flatnonzero(flagged).tolist() == BAD_PIXELS
    assert np.all(event.simulation.tel[1].true_image[BAD_PIXELS] == 0)
    assert event.simulation.tel[1].true_image.sum() == 1296 - len(BAD_PIXELS)
    # the R1 waveforms are kept, the charges of the bad pixels are interpolated
    np.testing.assert_array_equal(event.r1.tel[1].waveform, 1)
    # no bad pixels in tel 2
    assert not event.monitoring.tel[2].camera.coefficients.outlier_mask.any()
    assert event.simulation.tel[2].true_image.sum() == 1296


@pytest.mark.parametrize("calibrator_settings, expected", [
    ({}, [["type", "*", []], ["id", 22, [5, 7]]]),
    # the R0R1Calibrator section is used
    ({"R0R1Calibrator": {"bad_pixels": [["id", 21, [1]]]}}, [["id", 21, [1]]]),
])
def test_translate_legacy_bad_pixels(calibrator_settings, expected):

    legacy = {
        "analysis": {"bad_pixels": {"tel_021": [], "tel_022": [5, 7]}, "off_regions": 5},
        "telescope_calibration": {"tel_021": None, "tel_022": None, "bad_calib_px_interpolation": True},
        **calibrator_settings,
    }
    config = translate_legacy_calibration_config(legacy)

    assert config["analysis"] == {"off_regions": 5}
    assert config["R0R1Calibrator"]["bad_pixels"] == expected
    if not calibrator_settings:
        # the other legacy settings are translated too
        assert config["R0R1Calibrator"]["flag_bad_calibration_pixels"]


def test_translate_legacy_empty_bad_pixels():

    config = translate_legacy_calibration_config({"analysis": {"bad_pixels": {"tel_021": [], "tel_022": []}}})

    assert config == {"analysis": {}}


# window transmittance correction
@pytest.fixture
def window_file_of_ones(tmp_path):
    path = tmp_path / "window.txt"
    np.savetxt(path, np.column_stack([np.arange(1296), np.ones(1296)]), header="pixel_id\t correction_factor", comments="")
    return str(path)


def test_window_transmittance_correction(event, window_file_of_ones):

    waveforms = {}
    for name, settings in [("default", {}), ("ones", {"window_transmittance_file": window_file_of_ones})]:
        event.r1.tel.clear()
        calibrator(flag_bad_calibration_pixels=False, **settings)(without_coefficients(event), TEL_ID)
        waveforms[name] = event.r1.tel[TEL_ID].waveform.copy()

    factors = np.loadtxt((WINDOW_DIR / DEFAULT_WINDOW_FILES[TEL_ID]), skiprows=1, usecols=1)
    assert len(factors) == 1296 and not np.allclose(factors, 1)
    np.testing.assert_allclose(waveforms["default"], waveforms["ones"] / factors[:, np.newaxis])


def test_window_transmittance_file_per_telescope(window_file_of_ones):

    default = calibrator()
    assert default.window_transmittance_file_path(21).name == DEFAULT_WINDOW_FILES[21]
    assert default.window_transmittance_file_path(22).name == DEFAULT_WINDOW_FILES[22]

    custom = calibrator(window_transmittance_file=[["type", "*", None], ["id", 22, window_file_of_ones]])
    assert custom.window_transmittance_file_path(21).name == DEFAULT_WINDOW_FILES[21]
    np.testing.assert_array_equal(custom.window_transmittance(22), 1)


def test_window_transmittance_not_applied_to_simulations(mc_subarray):

    event = simulated_event(tel_ids=[1])
    mc_calibrator(mc_subarray)(event)

    np.testing.assert_array_equal(event.r1.tel[1].waveform, 1)


def test_translate_legacy_window_transmittance(window_file_of_ones):

    config = translate_legacy_calibration_config(
        {"window_transmittance": {"tel_021": None, "tel_022": window_file_of_ones}}
    )

    assert config == {"R0R1Calibrator": {"window_transmittance_file": [["type", "*", None], ["id", 22, window_file_of_ones]]}}
    assert translate_legacy_calibration_config({"window_transmittance": {"tel_021": None, "tel_022": None}}) == {}


# saturation correction
def old_saturated_charge_correction(event):
    """saturated_charge_correction before the ImageSaturationCorrector, as reference"""
    saturated_threshold, width_level, width_threshold, integration_level = 3000, 2500, 5, 0.2
    telescope = event.trigger.tels_with_trigger[0]
    r0data = event.r0.tel[telescope]
    waveforms = (r0data.waveform[0].T - r0data.pedestal)
    mask_saturated = np.max(waveforms, axis=0) > saturated_threshold
    saturated = False
    if sum(mask_saturated) > 0:
        image_new = event.dl1.tel[telescope].image
        peaktime_new = event.dl1.tel[telescope].peak_time
        for k, w in enumerate(waveforms.T):
            if mask_saturated[k]:
                mask_width = w > width_level
                n_samples = mask_width.shape[0]
                max_adc = max(w)
                mask_integration = w >= integration_level * max_adc
                integration_start = np.arange(0, n_samples)[mask_integration][0]
                index_max = np.arange(0, n_samples)[w == max_adc][0]
                int_stop = np.arange(0, n_samples)[(w < integration_level * max_adc) & (np.arange(0, n_samples) > index_max)]
                integration_stop = int_stop[0] if len(int_stop) > 0 else n_samples - 1
                width = sum(mask_width)
                peak_sample = width / 2 + np.arange(0, n_samples)[mask_width][0]
                if width > width_threshold:
                    image_new[k] = sum(event.r1.tel[telescope].waveform[0, k][integration_start:integration_stop + 1])
                    peaktime_new[k] = peak_sample * 4
                    saturated = True
        if saturated:
            event.dl1.tel[telescope].image = image_new
            event.dl1.tel[telescope].peak_time = peaktime_new
    return saturated


# saturated pixels: broad pulses (corrected), and one narrow pulse above the threshold (not corrected)
BROAD_PULSES = {10: (12, 22), 500: (5, 15), 1000: (35, 49)}
NARROW_PULSE = 700


def add_saturated_pulses(event, **calibrator_settings):
    """Event of the test file with saturated pulses added to the ADC samples, calibrated up to DL1"""
    event = deepcopy(event)
    waveform = event.r0.tel[TEL_ID].waveform
    pedestal = event.r0.tel[TEL_ID].pedestal
    for pixel, (start, stop) in BROAD_PULSES.items():
        waveform[0, pixel, start - 2:start] = pedestal[pixel] + 1200
        waveform[0, pixel, start:stop] = pedestal[pixel] + 3600
        waveform[0, pixel, stop:stop + 3] = pedestal[pixel] + 300
    waveform[0, NARROW_PULSE, 20:23] = pedestal[NARROW_PULSE] + 3600

    calibrator(**calibrator_settings)(without_coefficients(event), TEL_ID)
    CameraCalibrator(subarray=get_subarray(), config=CONFIG)(event)
    return event


@pytest.fixture
def saturated_event(event):
    return add_saturated_pulses(event)


def flagged_saturated(event):
    return np.flatnonzero(event.r1.tel[TEL_ID].pixel_status & PixelStatus.SATURATED).tolist()


def test_saturated_pixels_flagged_in_r1(saturated_event):

    # broad pulses, not the narrow one
    assert flagged_saturated(saturated_event) == sorted(BROAD_PULSES)
    assert np.all(saturated_event.r1.tel[TEL_ID].pixel_status & PixelStatus.HIGH_GAIN_STORED)


def test_image_saturation_corrector_as_before(saturated_event):

    reference = deepcopy(saturated_event)
    assert old_saturated_charge_correction(reference)

    corrector = ImageSaturationCorrector(subarray=get_subarray())
    image_before = saturated_event.dl1.tel[TEL_ID].image.copy()
    assert corrector(saturated_event, TEL_ID)

    dl1 = saturated_event.dl1.tel[TEL_ID]
    np.testing.assert_array_equal(dl1.image, reference.dl1.tel[TEL_ID].image)
    np.testing.assert_array_equal(dl1.peak_time, reference.dl1.tel[TEL_ID].peak_time)
    assert set(np.flatnonzero(dl1.image != image_before)) <= set(BROAD_PULSES)
    assert corrector.n_saturated_events == {TEL_ID: 1}


def test_image_saturation_corrector_peak_time(saturated_event):

    ImageSaturationCorrector(subarray=get_subarray())(saturated_event, TEL_ID)

    # middle of the samples above the width level, 4 ns per sample
    for pixel, (start, stop) in BROAD_PULSES.items():
        assert saturated_event.dl1.tel[TEL_ID].peak_time[pixel] == ((stop - start) / 2 + start) * 4


def test_saturation_settings(event):

    # all the pulses are narrower than 20 samples: no pixel flagged, nothing corrected
    saturated_event = add_saturated_pulses(event, saturation_width_threshold=20)
    image_before = saturated_event.dl1.tel[TEL_ID].image.copy()
    corrector = ImageSaturationCorrector(subarray=get_subarray())

    assert flagged_saturated(saturated_event) == []
    assert not corrector(saturated_event, TEL_ID)
    np.testing.assert_array_equal(saturated_event.dl1.tel[TEL_ID].image, image_before)
    assert corrector.n_saturated_events == {}


# calibration of each telescope set in the configuration (calibrate_r0, saturation levels)


def test_calibrate_r0_default_per_telescope(mc_subarray):
    """the R0 waveforms of the real telescopes are calibrated, not those of the simulated telescopes 1 and 2"""
    assert calibrator().calibrate_r0.tel[TEL_ID]
    assert not any(mc_calibrator(mc_subarray).calibrate_r0.tel[tel_id] for tel_id in (1, 2))


@pytest.mark.parametrize("nsb", ["low", "high"])
def test_no_saturation_correction_of_simulations(mc_subarray, nsb):
    """the MC configurations set no saturation correction (null levels) for the telescopes 1 and 2"""
    corrector = ImageSaturationCorrector(subarray=mc_subarray, config=load_config(MC_CONFIG_FILES[nsb], ismc=True))
    for tel_id in (1, 2):
        assert corrector.integration_level.tel[tel_id] is None
        assert corrector.peak_time_level.tel[tel_id] is None


def test_calibrate_r0_false_keeps_r1(event):
    """without R0 calibration, the R1 waveforms of the event are kept, the PDE drop is in the factor"""
    observed = without_coefficients(deepcopy(event))
    observed.r1.tel[TEL_ID] = R1CameraContainer(waveform=np.ones((1, 1296, 50), dtype=np.float32))
    calibrator(calibrate_r0=False, pde_drop_factor=0.5)(observed)
    np.testing.assert_array_equal(observed.r1.tel[TEL_ID].waveform, 1)
    np.testing.assert_allclose(observed.monitoring.tel[TEL_ID].camera.coefficients.factor, 2)


def test_calibrate_r0_of_simulation(mc_subarray):
    """the R0 waveforms of the simulations have no DigiCam pedestal: they cannot be calibrated"""
    event = simulated_event()
    event.r0.tel[1] = R0CameraContainer(waveform=np.ones((1, 1296, 50)))
    with pytest.raises(ValueError, match="calibrate_r0"):
        mc_calibrator(mc_subarray, calibrate_r0=True)(event)


@pytest.mark.parametrize("level", ["integration_level", "peak_time_level"])
def test_saturation_correction_disabled(saturated_event, level):
    """a null level: the saturated pixels are not corrected"""
    image_before = saturated_event.dl1.tel[TEL_ID].image.copy()
    corrector = ImageSaturationCorrector(
        subarray=get_subarray(), config=Config({"ImageSaturationCorrector": {level: None}})
    )
    assert not corrector(saturated_event, TEL_ID)
    np.testing.assert_array_equal(saturated_event.dl1.tel[TEL_ID].image, image_before)
    assert corrector.n_saturated_events == {}


def test_saturation_correction_with_factor(saturated_event):
    """the corrected charges are multiplied by the factor of the coefficients (voltage drop), as the others"""
    reference = deepcopy(saturated_event)
    ImageSaturationCorrector(subarray=get_subarray())(reference, TEL_ID)

    factor = np.linspace(1.1, 1.3, 1296)[np.newaxis]
    saturated_event.monitoring.tel[TEL_ID].camera.coefficients.factor = factor
    ImageSaturationCorrector(subarray=get_subarray())(saturated_event, TEL_ID)

    pixels = sorted(BROAD_PULSES)
    np.testing.assert_allclose(
        saturated_event.dl1.tel[TEL_ID].image[pixels], reference.dl1.tel[TEL_ID].image[pixels] * factor[0, pixels], rtol=1e-6,
    )
    np.testing.assert_array_equal(saturated_event.dl1.tel[TEL_ID].peak_time, reference.dl1.tel[TEL_ID].peak_time)


def test_saturated_outliers_not_corrected(saturated_event):
    """the saturated pixels which are outliers keep the charge interpolated by the CameraCalibrator"""
    outlier = sorted(BROAD_PULSES)[0]
    coefficients = saturated_event.monitoring.tel[TEL_ID].camera.coefficients
    coefficients.outlier_mask = coefficients.outlier_mask.copy()
    coefficients.outlier_mask[0, outlier] = True
    image_before = saturated_event.dl1.tel[TEL_ID].image.copy()
    ImageSaturationCorrector(subarray=get_subarray())(saturated_event, TEL_ID)

    corrected = np.flatnonzero(saturated_event.dl1.tel[TEL_ID].image != image_before)
    assert outlier not in corrected
    assert set(corrected) <= set(BROAD_PULSES) - {outlier}
