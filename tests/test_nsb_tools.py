import numpy as np
import pytest
import tables
from ctapipe.core import run_tool

from sst1mpipe.io import load_r0_pedestals
from sst1mpipe.resources import MONITORING_CONFIG_FILE, TEST_DATA_DIR
from sst1mpipe.scripts.sst1mpipe_monitoring_tool import MonitoringTool
from sst1mpipe.utils.NSB_tools import get_ped_table, mean_pedestal_std

# observation of Mrk 421 by tel 22: 13 pedestal events
INPUT_FILE = TEST_DATA_DIR / "zfits" / "SST1M2_20260120_1102.fits.fz"


@pytest.fixture(scope="module")
def monitoring_files(tmp_path_factory):
    """two monitoring files with chunks of 5 pedestal events, and a file without pedestal table"""
    directory = tmp_path_factory.mktemp("monitoring")
    paths = [directory / "run_1.monitoring.h5", directory / "run_2.monitoring.h5"]
    for path, max_events in zip(paths, [None, 20], strict=True):
        run_tool(MonitoringTool(), argv=[
            f"--input={INPUT_FILE}",
            f"--output={path}",
            f"--config={MONITORING_CONFIG_FILE}",
            "--no-progress",
            "--WaveformStatisticsCalculator.PedestalWaveformAggregator.SizeChunking.chunk_size=5",
            *([] if max_events is None else [f"--max-events={max_events}"]),
        ], raises=True)
    empty = directory / "empty.h5"
    tables.open_file(empty, "w").close()
    return paths, empty


@pytest.mark.parametrize("low_res", [False, True])
def test_get_ped_table(monitoring_files, low_res):
    paths, empty = monitoring_files
    pedestals = [load_r0_pedestals(path) for path in paths]
    assert len(pedestals[0]) > len(pedestals[1]) > 1

    # the files without pedestal table are skipped
    table = get_ped_table([empty, *paths], low_res=low_res)
    expected = [p[-1:] if low_res else p for p in pedestals]
    assert len(table) == sum(len(p) for p in expected)
    np.testing.assert_array_equal(table["event_id_end"], np.concatenate([p["event_id_end"] for p in expected]))
    assert np.all(np.isfinite(mean_pedestal_std(table)))

    assert get_ped_table([empty], low_res=low_res) is None
