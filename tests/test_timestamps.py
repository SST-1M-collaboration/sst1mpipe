import numpy as np
import pandas as pd
import pytest
from astropy.table import Table
from ctapipe.core import run_tool
from ctapipe.io import read_table

from sst1mpipe.io import add_time_ns
from sst1mpipe.io.sst1m_event_source import SST1MEventSource, camera_clock_to_time
from sst1mpipe.resources import RTA_CONFIG_FILE, TEST_DATA_DIR
from sst1mpipe.scripts.sst1mpipe_process_tool import ProcessorTool
from sst1mpipe.utils import get_timestamp_ns, time_to_ns

FILE_TEL_1 = (TEST_DATA_DIR / "zfits").joinpath('SST1M1_20260121_0001.fits.fz')
N_EVENTS = 20


def test_time_to_ns_camera_clock():
    # camera clocks (ns, TAI) one ns apart, not representable by a float64 unix time
    clock = np.uint64(1_768_953_600_123_456_789) + np.arange(5, dtype=np.uint64)
    time = camera_clock_to_time(clock)
    np.testing.assert_array_equal(time_to_ns(time), clock.astype(np.int64))
    assert time_to_ns(time[0]) == int(clock[0])


def test_get_timestamp_ns():
    times = np.array([1_768_953_600_123_456_789, 1_768_953_601_000_000_001], dtype=np.int64)
    assert np.array_equal(get_timestamp_ns(pd.DataFrame({'time_ns': times})), times)
    # DL1 files with the White Rabbit columns of older sst1mpipe versions
    legacy = Table({
        'time_wr_full_seconds': times // 1_000_000_000,
        'time_wr_frac_seconds': (times % 1_000_000_000) / 1e9,
    })
    assert np.all(np.abs(get_timestamp_ns(legacy) - times) <= 1)
    with pytest.raises(KeyError):
        get_timestamp_ns(pd.DataFrame({'local_time': [1.0]}))


@pytest.fixture(scope="module")
def dl1_file(tmp_path_factory):
    output = tmp_path_factory.mktemp("timestamps") / "events.dl1.h5"
    run_tool(ProcessorTool(), argv=[
        f"--input={FILE_TEL_1}",
        f"--output={output}",
        f"--config={RTA_CONFIG_FILE}",
        f"--max-events={N_EVENTS}",
        "--ProcessorTool.wobble_in_output_name=False",
    ], raises=True)
    return output


def test_add_time_ns(dl1_file):
    with SST1MEventSource([FILE_TEL_1], max_events=N_EVENTS) as source:
        clocks = {
            event.index.event_id: int(event.r0.tel[21].local_camera_clock)
            for event in source
        }

    parameters = read_table(dl1_file, "/dl1/event/telescope/parameters/tel_021")
    # the order of the rows is kept
    shuffled = parameters[np.random.default_rng(0).permutation(len(parameters))]
    params = add_time_ns(shuffled.copy(), str(dl1_file), 'tel_021')

    assert np.array_equal(params['event_id'], shuffled['event_id'])
    assert params['time_ns'].dtype == np.int64
    # trigger time of the telescope stored with the precision of the camera clock
    assert [clocks[e] for e in params['event_id']] == list(params['time_ns'])
    # added again, the column is replaced
    params = add_time_ns(params, str(dl1_file), 'tel_021')
    assert params.colnames.count('time_ns') == 1
