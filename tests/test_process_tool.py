import numpy as np
import pytest
import tables
from ctapipe.image.cleaning import NSBImageCleaner
from ctapipe.io import read_table
from ctapipe.io.hdf5dataformat import DL1_PIXEL_STATISTICS_GROUP

from conftest import OBSERVATION_RUNS, monitoring_source_args, run_process
from sst1mpipe.calib.calib import DEFAULT_CALIBRATION_FILES
from sst1mpipe.io import get_dl1_info, load_dl1_sst1m, pedestal_cleaning_fractions
from sst1mpipe.utils import get_survived_ped_fraction

TARGET, WOBBLE, RA, DEC = "MRK421", "W1", 166.9948, 38.1053


def process(run, output, *args, **kwargs):
    """sst1mpipe-process of the events of an observation run used in the tests"""
    return run_process(run["path"], output, f"--max-events={run['n_events']}", *args, **kwargs)


@pytest.mark.parametrize("tel_id", sorted(OBSERVATION_RUNS))
@pytest.mark.parametrize("flag_bad_calibration_pixels", [True, False])
def test_process_r0_file(tmp_path, observation_monitoring, tel_id, flag_bad_calibration_pixels):
    """regular processing of an observation run with its monitoring file, from R0 up to the DL1 parameters"""
    run = OBSERVATION_RUNS[tel_id]
    tool = process(
        run, tmp_path / "events.dl1.h5",
        *monitoring_source_args(observation_monitoring[tel_id]),
        f"--R0R1Calibrator.flag_bad_calibration_pixels={flag_bad_calibration_pixels}",
        wobble_in_output_name=True,
    )
    # the R0 -> R1 calibration is configured from the config file and the command line
    assert tool.r0_r1_calibrator.flag_bad_calibration_pixels.tel[tel_id] == flag_bad_calibration_pixels
    assert tool.n_skipped_events == 0
    assert tool.n_events_without_monitoring == 0

    # the wobble of the run is added to the name of the output file
    output = tmp_path / f"events_{WOBBLE}.dl1.h5"
    assert output.exists()

    # the pedestal events are not processed
    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
    assert len(parameters) == run["n_events"] - run["n_pedestals"]
    intensity = parameters["camera_frame_hillas_intensity"]
    assert np.isfinite(intensity).sum() > 0
    assert np.all(intensity[np.isfinite(intensity)] > 0)

    info = get_dl1_info(output)
    # calibration file of the telescope
    assert info["calib_file"][0].endswith(DEFAULT_CALIBRATION_FILES[tel_id])
    assert info[f"n_triggered_tel{tel_id - 20}"][0] == run["n_events"]
    assert info["n_pedestal"][0] == run["n_pedestals"]
    # the pedestal events are cleaned but not written: they do not survive the cleaning
    assert info["n_survived_pedestals"][0] == 0
    # fraction of pixels whose picture threshold is raised by the NSBImageCleaner
    assert 0 <= info["raised_threshold_fraction"][0] < 1
    assert pedestal_cleaning_fractions(output) == (0, info["raised_threshold_fraction"][0])
    zenith, _, ped_fractions, wobble, runs, _ = get_survived_ped_fraction([str(output)], tel=f"tel_{tel_id:03d}")
    assert list(ped_fractions) == [0] and list(wobble) == [WOBBLE]
    assert list(runs) == [int(run["path"].name.split(".")[0].split("_")[2])]
    assert 15 < zenith[0] < 17
    assert (info["target"][0], info["wobble"][0]) == (TARGET, WOBBLE)
    assert (info["ra"][0], info["dec"][0]) == pytest.approx((RA, DEC))
    assert not info["manual_coords"][0]


@pytest.fixture
def nsb_cleaner_pedestal_stds(monkeypatch):
    """the pedestal std given to the NSBImageCleaner for each event"""
    pedestal_stds = []
    nsb_image_cleaner_call = NSBImageCleaner.__call__

    def record(self, tel_id, image, arrival_times=None, *, monitoring=None):
        std = monitoring.pixel_statistics.pedestal_image.std
        pedestal_stds.append(None if std is None else np.array(std))
        return nsb_image_cleaner_call(self, tel_id, image, arrival_times, monitoring=monitoring)

    monkeypatch.setattr(NSBImageCleaner, "__call__", record)
    return pedestal_stds


def test_monitoring_used_by_nsb_image_cleaner(tmp_path, observation_monitoring, nsb_cleaner_pedestal_stds):
    tel_id, run = 21, OBSERVATION_RUNS[21]
    process(run, tmp_path / "events.dl1.h5", *monitoring_source_args(observation_monitoring[tel_id]))
    # std of the images of the single chunk of pedestal events (chunks of 1000 events), channel 0
    image = read_table(observation_monitoring[tel_id], f"{DL1_PIXEL_STATISTICS_GROUP}/sky_pedestal_image/tel_{tel_id:03d}")
    assert len(image) == 1
    # the images of all the events, pedestal events included
    assert len(nsb_cleaner_pedestal_stds) == run["n_events"]
    for std in nsb_cleaner_pedestal_stds:
        np.testing.assert_allclose(std, image["std"][0, 0])
    assert 0.1 < np.median(image["std"][0, 0]) < 2


def test_process_without_monitoring_source(tmp_path, nsb_cleaner_pedestal_stds):
    """the raw data are calibrated without the statistics of the pedestal events, and without pointing"""
    tel_id, run = 22, OBSERVATION_RUNS[22]
    output, log_file = tmp_path / "events.dl1.h5", tmp_path / "process.log"
    tool = process(run, output, f"--log-file={log_file}")
    assert "Telescope 22: no camera calibration coefficients" in log_file.read_text()
    assert tool.monitoring_sources == []

    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
    assert len(parameters) == run["n_events"] - run["n_pedestals"]
    assert np.isfinite(parameters["camera_frame_hillas_intensity"]).sum() > 0
    info = get_dl1_info(output)
    assert info["n_pedestal"][0] == run["n_pedestals"]
    # no statistics of the pedestal events for the NSBImageCleaner: no raised picture threshold
    assert nsb_cleaner_pedestal_stds == [None] * run["n_events"]
    assert np.isnan(info["raised_threshold_fraction"][0])

    with tables.open_file(output) as dl1:
        assert "/dl0/monitoring" not in dl1 and "/r0/monitoring" not in dl1 and "/dl1/monitoring" not in dl1
    # without pointing, the DL1 file cannot be loaded for the reconstruction (DL2)
    with pytest.raises(SystemExit):
        load_dl1_sst1m(str(output), tel=f"tel_{tel_id:03d}")
