import numpy as np
import pytest
import tables

from conftest import DARK_FILE, OBSERVATION_RUNS, WAVEFORM_CHUNK_ARGS, run_monitoring
from sst1mpipe.io import load_r0_pedestals
from sst1mpipe.io.containers import PEDESTAL_EVENT_TYPES
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils.NSB_tools import get_dark_baseline, get_ped_table, mean_pedestal_std


@pytest.fixture(scope="module")
def monitoring_files(observation_monitoring, tmp_path_factory):
    """
    monitoring files of the observation of tel 22 and of its 20 first events, with the chunks of
    WAVEFORM_CHUNK_ARGS, and a file without pedestal table
    """
    directory = tmp_path_factory.mktemp("monitoring")
    first_events = directory / "first_events.monitoring.h5"
    run_monitoring(
        OBSERVATION_RUNS[22]["path"], first_events, "--max-events=20",
        *WAVEFORM_CHUNK_ARGS,
    )
    empty = directory / "empty.h5"
    tables.open_file(empty, "w").close()
    return [observation_monitoring[22], first_events], empty


@pytest.mark.parametrize("low_res", [False, True])
def test_get_ped_table(monitoring_files, low_res):
    paths, empty = monitoring_files
    pedestals = [load_r0_pedestals(path) for path in paths]
    assert len(pedestals[0]) > len(pedestals[1]) > 1

    # the tables of the files sorted by name; the files without pedestal table are skipped
    table = get_ped_table([empty, *paths], low_res=low_res)
    expected = [p[-1:] if low_res else p for _, p in sorted(zip(map(str, paths), pedestals, strict=True))]
    assert len(table) == sum(len(p) for p in expected)
    np.testing.assert_array_equal(table["event_id_end"], np.concatenate([p["event_id_end"] for p in expected]))
    assert np.all(np.isfinite(mean_pedestal_std(table)))

    assert get_ped_table([empty], low_res=low_res) is None


def test_get_dark_baseline():
    """mean of the first 50 ADC samples of each pixel over the pedestal events of the dark run"""
    max_evt = 100
    with SST1MEventSource(DARK_FILE, max_events=max_evt) as source:
        samples = np.array([
            event.r0.tel[event.trigger.tels_with_trigger[0]].waveform[0][:, :50]
            for event in source if event.trigger.event_type in PEDESTAL_EVENT_TYPES
        ], dtype=np.float64)
    assert len(samples) > 10
    baseline = get_dark_baseline(DARK_FILE, max_evt=max_evt)
    assert baseline.shape == (1296, )
    np.testing.assert_allclose(baseline, samples.mean(axis=(0, 2)), rtol=1e-12)
    # no event of the requested types
    with pytest.raises(ValueError, match="No event"):
        get_dark_baseline(DARK_FILE, max_evt=max_evt, event_types=())
