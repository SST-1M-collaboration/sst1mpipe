"""
Shared test data: the SST-1M runs of the test data, the tools run on them (sst1mpipe-monitoring,
sst1mpipe-process) and a reader of the events of the zfits files.
The constants and functions are imported by the test modules (``from conftest import ...``).
"""
import numpy as np
import pytest
from astropy.utils import iers
from ctapipe.core import run_tool
from protozfits import File

from sst1mpipe.resources import MONITORING_CONFIG_FILE, RTA_CONFIG_FILE, TEST_DATA_DIR
from sst1mpipe.scripts.sst1mpipe_monitoring_tool import MonitoringTool
from sst1mpipe.scripts.sst1mpipe_process_tool import ProcessorTool

# the tests do not download the IERS tables (Earth orientation): the alt/az of the pointing
# are computed with the tables of astropy, precise enough for the tests (< 1 arcmin), without
# depending on the network
iers.conf.auto_download = False

ZFITS_DIR = TEST_DATA_DIR / "zfits"
S_TO_NS = 1_000_000_000

# observation of Mrk 421 (wobble W1) by both telescopes, with Cherenkov and pedestal events:
# the run of tel 22 is taken during the run of tel 21. Number of events (and of pedestal events)
# used in the tests: the 150 first events of tel 21, all the events of tel 22
OBSERVATION_RUNS = {
    21: dict(path=ZFITS_DIR / "SST1M1_20260120_1179.fits.fz", n_events=150, n_pedestals=40),
    22: dict(path=ZFITS_DIR / "SST1M2_20260120_1102.fits.fz", n_events=35, n_pedestals=13),
}
# dark run of tel 22: pedestal events only
DARK_FILE = ZFITS_DIR / "SST1M2_20260119_0007.fits.fz"
# runs of the transition between two wobbles (night of 2026-01-21)
TRANSITION_FILES = {
    21: ZFITS_DIR / "SST1M1_20260121_0206.fits.fz",
    22: ZFITS_DIR / "SST1M2_20260121_0585.fits.fz",
}

# chunks of pedestal events of the statistics of the ADC samples, in the monitoring files of the tests
WAVEFORM_CHUNK_SIZE = 5


def run_monitoring(input_file, output, *args):
    """sst1mpipe-monitoring with the monitoring config, returns the tool"""
    tool = MonitoringTool()
    run_tool(tool, argv=[
        f"--input={input_file}",
        f"--output={output}",
        f"--config={MONITORING_CONFIG_FILE}",
        "--no-progress",
        *args,
    ], raises=True)
    return tool


def run_process(input_file, output, *args, config=RTA_CONFIG_FILE, wobble_in_output_name=False):
    """
    sst1mpipe-process (RTA config by default), returns the tool.
    The output file keeps its name, unless wobble_in_output_name.
    """
    tool = ProcessorTool()
    run_tool(tool, argv=[
        f"--input={input_file}",
        f"--output={output}",
        f"--config={config}",
        "--no-progress",
        f"--ProcessorTool.wobble_in_output_name={wobble_in_output_name}",
        *args,
    ], raises=True)
    return tool


def monitoring_source_args(monitoring_file):
    """sst1mpipe-process arguments to use a monitoring file of sst1mpipe-monitoring"""
    return ["--monitoring-source=SST1MMonitoringSource", f"--SST1MMonitoringSource.input_files={monitoring_file}"]


def read_zfits_events(path, events=slice(None)):
    """
    SWAT event id, White Rabbit time (s, ns and total ns) and camera event type of the events of
    a zfits file, read with protozfits
    """
    with File(str(path)) as f:
        # the slice is limited to the events of the file
        rows = list(f.Events[slice(*events.indices(len(f.Events)))])
    time_s = np.array([row.local_time_sec for row in rows], dtype=np.int64)
    time_ns = np.array([row.local_time_nanosec for row in rows], dtype=np.int64)
    return dict(
        event_id=np.array([row.arrayEvtNum for row in rows], dtype=np.int64),
        time_s=time_s,
        time_ns=time_ns,
        time=time_s * S_TO_NS + time_ns,
        event_type=np.array([row.event_type for row in rows]),
    )


@pytest.fixture(scope="session")
def observation_monitoring(tmp_path_factory):
    """
    monitoring files of the observation runs (OBSERVATION_RUNS) by telescope id,
    with chunks of WAVEFORM_CHUNK_SIZE pedestal events
    """
    directory = tmp_path_factory.mktemp("observation_monitoring")
    paths = {}
    for tel_id, run in OBSERVATION_RUNS.items():
        paths[tel_id] = directory / f"{run['path'].name.split('.')[0]}.monitoring.h5"
        run_monitoring(
            run["path"], paths[tel_id],
            f"--max-events={run['n_events']}",
            f"--WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_size={WAVEFORM_CHUNK_SIZE}",
        )
    return paths
