import numpy as np
import pytest
import tables
from ctapipe.core import run_tool
from ctapipe.image.cleaning import NSBImageCleaner
from ctapipe.io import read_table
from ctapipe.io.hdf5dataformat import DL1_PIXEL_STATISTICS_GROUP

from sst1mpipe.calib.calib import DEFAULT_CALIBRATION_FILES
from sst1mpipe.io import get_dl1_info, load_dl1_sst1m
from sst1mpipe.scripts.sst1mpipe_monitoring_tool import MonitoringTool
from sst1mpipe.scripts.sst1mpipe_process_tool import ProcessorTool
from sst1mpipe.resources import MONITORING_CONFIG_FILE, RTA_CONFIG_FILE, TEST_DATA_DIR

# observation of Mrk 421 (wobble W1) by both telescopes, with Cherenkov and pedestal events
OBSERVATION_FILES = {
    21: dict(path=TEST_DATA_DIR / "zfits" / "SST1M1_20260120_1179.fits.fz", n_events=150, n_pedestals=40),
    22: dict(path=TEST_DATA_DIR / "zfits" / "SST1M2_20260120_1102.fits.fz", n_events=35, n_pedestals=13),
}
TARGET, WOBBLE, RA, DEC = "MRK421", "W1", 166.9948, 38.1053


def process(run, output, *args):
    tool = ProcessorTool()
    run_tool(tool, argv=[
        f"--input={run['path']}",
        f"--output={output}",
        f"--config={RTA_CONFIG_FILE}",
        f"--max-events={run['n_events']}",
        "--no-progress",
        *args,
    ], raises=True)
    return tool


@pytest.fixture(scope="module")
def monitoring_files(tmp_path_factory):
    """monitoring files of the first events of the runs (sst1mpipe-monitoring), with chunks of 10 pedestal events"""
    directory = tmp_path_factory.mktemp("monitoring")
    paths = {}
    for tel_id, run in OBSERVATION_FILES.items():
        paths[tel_id] = directory / f"tel_{tel_id}.monitoring.h5"
        run_tool(MonitoringTool(), argv=[
            f"--input={run['path']}",
            f"--output={paths[tel_id]}",
            f"--config={MONITORING_CONFIG_FILE}",
            f"--max-events={run['n_events']}",
            "--no-progress",
            "--WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_size=10",
        ], raises=True)
    return paths


def monitoring_source(monitoring_file):
    return ["--monitoring-source=SST1MMonitoringSource", f"--SST1MMonitoringSource.input_files={monitoring_file}"]


@pytest.mark.parametrize("tel_id", sorted(OBSERVATION_FILES))
@pytest.mark.parametrize("voltage_drop_correction", ["global", "none"])
def test_process_r0_file(tmp_path, monitoring_files, tel_id, voltage_drop_correction):
    """regular processing of an observation run with its monitoring file, from R0 up to the DL1 parameters"""
    run = OBSERVATION_FILES[tel_id]
    tool = process(
        run, tmp_path / "events.dl1.h5",
        *monitoring_source(monitoring_files[tel_id]),
        f"--R0R1Calibrator.voltage_drop_correction={voltage_drop_correction}",
    )
    # the R0 -> R1 calibration is configured from the config file and the command line
    assert tool.r0_r1_calibrator.voltage_drop_correction.tel[tel_id] == voltage_drop_correction
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
    # the pedestal events are not processed
    assert np.isnan(info["n_survived_pedestals"][0])
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


def test_monitoring_used_by_nsb_image_cleaner(tmp_path, monitoring_files, nsb_cleaner_pedestal_stds):
    tel_id, run = 21, OBSERVATION_FILES[21]
    process(
        run, tmp_path / "events.dl1.h5",
        *monitoring_source(monitoring_files[tel_id]),
        "--ProcessorTool.wobble_in_output_name=False",
    )
    # std of the images of the single chunk of pedestal events (chunks of 1000 events), channel 0
    image = read_table(monitoring_files[tel_id], f"{DL1_PIXEL_STATISTICS_GROUP}/sky_pedestal_image/tel_{tel_id:03d}")
    assert len(image) == 1
    assert len(nsb_cleaner_pedestal_stds) == run["n_events"] - run["n_pedestals"]
    for std in nsb_cleaner_pedestal_stds:
        np.testing.assert_allclose(std, image["std"][0, 0])
    assert 0.1 < np.median(image["std"][0, 0]) < 2


def test_process_without_monitoring_source(tmp_path, nsb_cleaner_pedestal_stds):
    """the raw data are calibrated without the statistics of the pedestal events, and without pointing"""
    tel_id, run = 22, OBSERVATION_FILES[22]
    output, log_file = tmp_path / "events.dl1.h5", tmp_path / "process.log"
    tool = process(run, output, "--ProcessorTool.wobble_in_output_name=False", f"--log-file={log_file}")
    assert "No monitoring source" in log_file.read_text()
    assert tool.monitoring_sources == []

    parameters = read_table(output, f"/dl1/event/telescope/parameters/tel_{tel_id:03d}")
    assert len(parameters) == run["n_events"] - run["n_pedestals"]
    assert np.isfinite(parameters["camera_frame_hillas_intensity"]).sum() > 0
    assert get_dl1_info(output)["n_pedestal"][0] == run["n_pedestals"]
    # no statistics of the pedestal events for the NSBImageCleaner
    assert nsb_cleaner_pedestal_stds == [None] * len(parameters)

    with tables.open_file(output) as dl1:
        assert "/dl0/monitoring" not in dl1 and "/r0/monitoring" not in dl1 and "/dl1/monitoring" not in dl1
    # without pointing, the DL1 file cannot be loaded for the reconstruction (DL2)
    with pytest.raises(SystemExit):
        load_dl1_sst1m(str(output), tel=f"tel_{tel_id:03d}")
