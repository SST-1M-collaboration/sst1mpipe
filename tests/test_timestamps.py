import shutil
from itertools import islice

import numpy as np
import pytest
import tables
from ctapipe.core import run_tool
from ctapipe.io import read_table
from ctapipe.time import time_to_ctao_high_res
from protozfits import File

from sst1mpipe.io import load_dl1_sst1m, read_trigger_time
from sst1mpipe.resources import RTA_CONFIG_FILE, TEST_DATA_DIR
from sst1mpipe.scripts.sst1mpipe_process_tool import ProcessorTool
from sst1mpipe.utils import get_trigger_time_ns

N_EVENTS = 50
S_TO_NS = 1_000_000_000
FILES = {
    21: (TEST_DATA_DIR / "zfits").joinpath('SST1M1_20260121_0001.fits.fz'),
    22: (TEST_DATA_DIR / "zfits").joinpath('SST1M2_20260121_0001.fits.fz'),
}


def read_zfits_times(path, n_events):
    """event id (SWAT) and White Rabbit time (s, ns) of the first events of the zfits file"""
    with File(str(path)) as f:
        events = list(islice(f.Events, n_events))
    return (
        np.array([e.arrayEvtNum for e in events], dtype=np.int64),
        np.array([e.local_time_sec for e in events], dtype=np.int64),
        np.array([e.local_time_nanosec for e in events], dtype=np.int64),
    )


@pytest.fixture(scope="module", params=sorted(FILES))
def dl1_file(request, tmp_path_factory):
    tel_id = request.param
    output = tmp_path_factory.mktemp(f"timestamps_{tel_id}") / "events.dl1.h5"
    run_tool(ProcessorTool(), argv=[
        f"--input={FILES[tel_id]}",
        f"--output={output}",
        f"--config={RTA_CONFIG_FILE}",
        f"--max-events={N_EVENTS}",
        "--ProcessorTool.wobble_in_output_name=False",
        # the test files are dark runs
        "--ProcessorTool.allowed_sb_types=CALIBRATION",
    ], raises=True)
    return tel_id, output


@pytest.mark.parametrize("table", ["/dl1/event/subarray/trigger", "/dl1/event/telescope/trigger"])
def test_white_rabbit_time_stored_in_dl1(dl1_file, table):
    """the stored time is exactly the White Rabbit time of the zfits file: seconds and ns"""
    tel_id, output = dl1_file
    event_id, seconds, nanoseconds = read_zfits_times(FILES[tel_id], N_EVENTS)

    with tables.open_file(output) as h5:
        rows = h5.get_node(table).read()
        assert h5.get_node(table).attrs["CTAFIELD_2_TIME_SCALE" if "subarray" in table
                                        else "CTAFIELD_3_TIME_SCALE"] == "tai"

    np.testing.assert_array_equal(rows["event_id"], event_id)
    # ctao_high_res format: (seconds, quarter of nanoseconds)
    np.testing.assert_array_equal(rows["time"][:, 0], seconds)
    np.testing.assert_array_equal(rows["time"][:, 1], 4 * nanoseconds)


def test_white_rabbit_time_read_from_dl1(dl1_file):
    """the time read by ctapipe from the DL1 file keeps the ns precision"""
    tel_id, output = dl1_file
    _, seconds, nanoseconds = read_zfits_times(FILES[tel_id], N_EVENTS)
    time_ns = seconds * S_TO_NS + nanoseconds

    time = read_table(output, "/dl1/event/subarray/trigger")["time"]
    assert time.scale == "tai"

    # astropy keeps the time as two floats: the conversion of a time difference to a single
    # float is precise to ~0.01 ns, so the times are compared rounded to the ns
    from_start = (time - time[0]).to_value("ns")
    np.testing.assert_array_equal(np.round(from_start).astype(np.int64), time_ns - time_ns[0])

    # time differences between the events, to the ns
    np.testing.assert_array_equal(
        np.round((time[1:] - time[:-1]).to_value("ns")).astype(np.int64), np.diff(time_ns)
    )
    # absolute times: converted back to integers (seconds, quarter of ns) by ctapipe
    high_res = time_to_ctao_high_res(time).astype(np.int64)
    np.testing.assert_array_equal(high_res[:, 0] * S_TO_NS + high_res[:, 1] // 4, time_ns)


def test_read_trigger_time(dl1_file):
    """the trigger times read from the DL1 file are the White Rabbit times of the zfits file"""
    tel_id, output = dl1_file
    event_id, seconds, nanoseconds = read_zfits_times(FILES[tel_id], N_EVENTS)

    trigger_time = read_trigger_time(output, f"tel_{tel_id:03d}")
    assert trigger_time["trigger_time"].scale == "tai"
    np.testing.assert_array_equal(trigger_time["event_id"], event_id)
    np.testing.assert_array_equal(get_trigger_time_ns(trigger_time), seconds * S_TO_NS + nanoseconds)
    # no event of the other telescope
    assert len(read_trigger_time(output, "tel_021" if tel_id == 22 else "tel_022")) == 0


@pytest.mark.parametrize("table", ["astropy", "pandas"])
def test_load_dl1_sst1m_trigger_time(dl1_file, tmp_path, table):
    """load_dl1_sst1m adds the trigger time of each event to the DL1 table, with ns precision"""
    tel_id, output = dl1_file
    tel = f"tel_{tel_id:03d}"
    event_id, seconds, nanoseconds = read_zfits_times(FILES[tel_id], N_EVENTS)
    clock = dict(zip(event_id, seconds * S_TO_NS + nanoseconds, strict=True))

    # parameters table in another order than the trigger table, with the pointing
    # columns written by write_extra_parameters
    dl1 = tmp_path / "events.dl1.h5"
    shutil.copy(output, dl1)
    params = read_table(dl1, f"/dl1/event/telescope/parameters/{tel}")
    params = params[np.random.default_rng(0).permutation(len(params))]
    params["true_az_tel"] = 0.0
    params["true_alt_tel"] = 70.0
    params.write(dl1, path=f"/dl1/event/telescope/parameters/{tel}", overwrite=True, append=True)

    data = load_dl1_sst1m(str(dl1), tel=tel, table=table)
    np.testing.assert_array_equal(data["event_id"], params["event_id"])
    assert list(get_trigger_time_ns(data)) == [clock[e] for e in data["event_id"]]
