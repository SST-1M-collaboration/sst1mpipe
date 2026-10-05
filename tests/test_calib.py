from importlib.resources import files

import numpy as np
import pytest
import json

from ctapipe.calib import CameraCalibrator
from ctapipe.containers import ArrayEventContainer, R1CameraContainer, SimulatedEventContainer
from ctapipe.instrument import SubarrayDescription
from traitlets.config import Config

from sst1mpipe.calib import R0R1Calibrator, saturated_charge_correction
from sst1mpipe.calib.calib import DEFAULT_CALIBRATION_FILES
from sst1mpipe.io import load_config, translate_legacy_calibration_config
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils import get_subarray

FILE_TEL_1 = files('sst1mpipe.resources.zfits').joinpath('SST1M1_20260121_0001.fits.fz')
DATA_CONFIG_FILE = files('sst1mpipe.data').joinpath('sst1mpipe_data_config.json')
MC_CONFIG_FILE = files('sst1mpipe.data').joinpath('sst1mpipe_mc_config.json')
CONFIG = load_config(DATA_CONFIG_FILE, ismc=False)
TEL_ID = 21
CALIBRATION_FILE_TEL_2 = str(files('sst1mpipe.data').joinpath(DEFAULT_CALIBRATION_FILES[22]))


@pytest.fixture(scope="module")
def event():
    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=1)
    return next(iter(source))


def calibrator(**settings):
    return R0R1Calibrator(subarray=get_subarray(), config=Config({"R0R1Calibrator": settings}))


def with_pedestal_std(event, std):
    event.mon.tel[TEL_ID].r0.charge_std = std
    return event


def test_r0_r1_dl1_calibration():

    source = SST1MEventSource(input_url=FILE_TEL_1, max_events=3)
    n_pixels = source.subarray.tel[TEL_ID].camera.readout.n_pixels
    calibrator_r0_r1 = R0R1Calibrator(subarray=source.subarray, config=CONFIG)
    r1_dl1_calibrator = CameraCalibrator(subarray=source.subarray, config=CONFIG)

    for event in source:
        r0_waveform = event.r0.tel[TEL_ID].waveform.copy()

        calibrator_r0_r1(event, TEL_ID)
        r1 = event.r1.tel[TEL_ID]

        # ctapipe shape: (n_channels, n_pixels, n_samples)
        assert r1.waveform.shape == r0_waveform.shape == (1, n_pixels, r0_waveform.shape[-1])
        bad_pixels = event.mon.tel[TEL_ID].pixel_status.hardware_failing_pixels[0]
        assert bad_pixels.any()
        assert calibrator_r0_r1.n_bad_pixels[TEL_ID] == bad_pixels.sum()
        assert np.all(r1.waveform[:, bad_pixels] == 0)
        # the raw data are not modified by the calibration
        np.testing.assert_array_equal(event.r0.tel[TEL_ID].waveform, r0_waveform)

        r1_dl1_calibrator(event)
        assert event.dl1.tel[TEL_ID].image.shape == (n_pixels, )
        assert isinstance(saturated_charge_correction(event), bool)


def test_all_telescopes_with_r0_data_are_calibrated(event):

    event.r1.tel.clear()
    calibrator()(event)

    assert list(event.r1.tel.keys()) == [TEL_ID]


def test_pedestal_subtraction_and_dc_to_pe(event):

    calibrator_r0_r1 = calibrator(voltage_drop_correction="none", flag_bad_calibration_pixels=False)
    calibrator_r0_r1(with_pedestal_std(event, None), TEL_ID)

    r0 = event.r0.tel[TEL_ID]
    dc_to_pe, _ = calibrator_r0_r1.calibration_parameters(TEL_ID)
    expected = (r0.waveform[0] - r0.pedestal[:, np.newaxis]) / dc_to_pe[:, np.newaxis]
    np.testing.assert_allclose(event.r1.tel[TEL_ID].waveform[0], expected)
    assert calibrator_r0_r1.n_bad_pixels[TEL_ID] == 0


@pytest.mark.parametrize("correction", ["none", "global", "pixelwise"])
def test_voltage_drop_correction(event, correction):

    pedestal_std = np.linspace(3, 6, 1296)
    no_correction = calibrator(voltage_drop_correction="none", flag_dead_pixels=False)
    no_correction(with_pedestal_std(event, pedestal_std), TEL_ID)
    waveform = event.r1.tel[TEL_ID].waveform.copy()

    calibrator_r0_r1 = calibrator(voltage_drop_correction=correction, flag_dead_pixels=False)
    calibrator_r0_r1(with_pedestal_std(event, pedestal_std), TEL_ID)

    voltage_drop = np.broadcast_to(calibrator_r0_r1.voltage_drop(TEL_ID, pedestal_std), (1296, ))
    if correction == "none":
        assert np.all(voltage_drop == 1)
    elif correction == "global":
        # charges divided by the drop of the current, smaller for a higher NSB (pedestal variance)
        assert np.all(voltage_drop == voltage_drop[0]) and voltage_drop[0] < 1
    else:
        assert np.all(np.diff(voltage_drop) < 0)
    np.testing.assert_allclose(event.r1.tel[TEL_ID].waveform[0], waveform[0] / voltage_drop[:, np.newaxis])


def test_no_voltage_drop_correction_without_pedestals(event):

    calibrator_r0_r1 = calibrator(voltage_drop_correction="pixelwise")

    assert calibrator_r0_r1.voltage_drop(TEL_ID, None) == 1.0


@pytest.mark.parametrize("flag_bad_calibration, flag_dead, threshold, expected_dead", [
    (False, False, 2.5, False),
    (True, False, 2.5, False),
    (False, True, 2.5, True),
    (True, True, 2.5, True),
    (True, True, 0.5, False),
])
def test_bad_pixels(event, flag_bad_calibration, flag_dead, threshold, expected_dead):

    dead_pixel = 100
    pedestal_std = np.full((1, 1296), 4.0)
    pedestal_std[0, dead_pixel] = 1.0
    calibrator_r0_r1 = calibrator(
        flag_bad_calibration_pixels=flag_bad_calibration,
        flag_dead_pixels=flag_dead,
        dead_pixel_std_threshold=threshold,
    )
    calibrator_r0_r1(with_pedestal_std(event, pedestal_std), TEL_ID)

    _, mask_bad_calibration = calibrator_r0_r1.calibration_parameters(TEL_ID)
    flagged = event.mon.tel[TEL_ID].pixel_status.hardware_failing_pixels[0]
    assert mask_bad_calibration.any() and not mask_bad_calibration[dead_pixel]
    assert np.all(flagged[mask_bad_calibration] == flag_bad_calibration)
    assert flagged[dead_pixel] == expected_dead
    assert np.all(event.r1.tel[TEL_ID].waveform[:, flagged] == 0)
    assert np.all(event.r1.tel[TEL_ID].waveform[:, ~flagged].any(axis=-1))


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


@pytest.mark.parametrize("settings", [
    {"voltage_drop_correction": "wrong"},
    {"calibration_file": "/does/not/exist.h5"},
])
def test_wrong_settings(settings):

    with pytest.raises(Exception, match="(?i)trait|exist"):
        calibrator(**settings)


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
        "flag_dead_pixels": False,
        "voltage_drop_correction": expected,
        "mc_pde_correction": False,
    }
    calibrator_r0_r1 = R0R1Calibrator(subarray=get_subarray(), config=Config(config))
    assert calibrator_r0_r1.voltage_drop_correction.tel[21] == expected
    assert str(calibrator_r0_r1.calibration_file_path(22)) == CALIBRATION_FILE_TEL_2


def test_config_without_legacy_settings_is_unchanged():

    config = {"R0R1Calibrator": {"voltage_drop_correction": "none"}, "NsbCalibrator": {}}

    assert translate_legacy_calibration_config(config) == config


def test_default_config_settings():

    calibrator_r0_r1 = R0R1Calibrator(subarray=get_subarray(), config=CONFIG)

    for tel_id in (21, 22):
        assert calibrator_r0_r1.voltage_drop_correction.tel[tel_id] == "global"
        assert calibrator_r0_r1.flag_bad_calibration_pixels.tel[tel_id]
        assert calibrator_r0_r1.flag_dead_pixels.tel[tel_id]
        assert calibrator_r0_r1.dead_pixel_std_threshold.tel[tel_id] == 2.5
        assert calibrator_r0_r1.calibration_file_path(tel_id).name == DEFAULT_CALIBRATION_FILES[tel_id]


# PDE drop correction of the simulations
PDE_FILES = ["qe_SST1M_5477_ave_TEL1_NSB251.0", "qe_SST1M_5477_ave_TEL2_NSB300.0", "qe_dummy"]
PDE_DROP_FACTORS = {1: 0.9393873691700036, 2: 0.9819055121768047}  # mc_pde_correction_factors.json


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


def mc_calibrator(mc_subarray, simulated_pde_files=PDE_FILES, **settings):
    return R0R1Calibrator(
        subarray=mc_subarray, config=Config({"R0R1Calibrator": settings}),
        simulated_pde_files=simulated_pde_files,
    )


def test_mc_pde_correction(mc_subarray):

    event = simulated_event()
    calibrator_r0_r1 = mc_calibrator(mc_subarray)
    calibrator_r0_r1(event)

    for tel_id, factor in PDE_DROP_FACTORS.items():
        assert calibrator_r0_r1.pde_drop_factor(tel_id) == factor
        np.testing.assert_allclose(event.r1.tel[tel_id].waveform, 1 / factor, rtol=1e-6)
    # the R0 data of the simulations are not used
    assert len(event.r0.tel) == 0


def test_mc_pde_correction_single_telescope(mc_subarray):

    event = simulated_event()
    mc_calibrator(mc_subarray)(event, 2)

    np.testing.assert_array_equal(event.r1.tel[1].waveform, 1)
    np.testing.assert_allclose(event.r1.tel[2].waveform, 1 / PDE_DROP_FACTORS[2], rtol=1e-6)


def test_mc_pde_correction_disabled(mc_subarray):

    event = simulated_event()
    # PDE files are not needed if the correction is disabled
    calibrator_r0_r1 = mc_calibrator(
        mc_subarray, simulated_pde_files=None, mc_pde_correction=[["type", "*", True], ["id", 1, False]],
    )
    calibrator_r0_r1(event, 1)
    np.testing.assert_array_equal(event.r1.tel[1].waveform, 1)

    with pytest.raises(ValueError, match="simulated_pde_files"):
        calibrator_r0_r1(event, 2)


def test_mc_pde_correction_unknown_pde_file(mc_subarray):

    calibrator_r0_r1 = mc_calibrator(mc_subarray, simulated_pde_files=["qe_unknown"])

    with pytest.raises(ValueError, match="No PDE drop correction factor of telescope 1"):
        calibrator_r0_r1(simulated_event(tel_ids=[1]))


def test_mc_pde_correction_file(mc_subarray, tmp_path):

    path = tmp_path / "pde_factors.json"
    path.write_text(json.dumps({"mc_correction_for_PDE": {"tel_001": {"qe_custom": 0.5}}}))
    event = simulated_event(tel_ids=[1])
    mc_calibrator(mc_subarray, simulated_pde_files=["qe_custom"], mc_pde_correction_file=str(path))(event)

    np.testing.assert_array_equal(event.r1.tel[1].waveform, 2)


def test_mc_config_settings(mc_subarray):

    config = load_config(MC_CONFIG_FILE, ismc=True)
    calibrator_r0_r1 = R0R1Calibrator(subarray=mc_subarray, config=config, simulated_pde_files=PDE_FILES)

    assert "mc_correction_for_PDE" not in config["NsbCalibrator"]
    assert "intensity_correction" in config["NsbCalibrator"]  # used by sst1mpipe_dl1_dl2
    assert calibrator_r0_r1.mc_pde_correction.tel[1] and calibrator_r0_r1.mc_pde_correction.tel[2]
    assert calibrator_r0_r1.pde_drop_factor(1) == PDE_DROP_FACTORS[1]


def test_data_events_are_not_pde_corrected(event):

    # observed data: R0 -> R1 calibration, no PDE drop correction (the PDE files are not needed)
    event.r1.tel.clear()
    calibrator(voltage_drop_correction="none", flag_bad_calibration_pixels=False)(
        with_pedestal_std(event, None), TEL_ID,
    )
    r0 = event.r0.tel[TEL_ID]
    dc_to_pe, _ = calibrator().calibration_parameters(TEL_ID)
    np.testing.assert_allclose(
        event.r1.tel[TEL_ID].waveform[0], (r0.waveform[0] - r0.pedestal[:, np.newaxis]) / dc_to_pe[:, np.newaxis],
    )
