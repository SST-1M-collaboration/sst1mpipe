"""sst1mpipe-monitoring: monitoring data of the SST-1M runs"""
import astropy.units as u
import numpy as np
import pytest
import tables
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time
from ctapipe.containers import ArrayEventContainer, SchedulingBlockType
from ctapipe.core.tool import ToolConfigurationError
from ctapipe.instrument import SubarrayDescription
from ctapipe.io import HDF5MonitoringSource, read_table
from ctapipe.io.hdf5dataformat import DL0_TEL_POINTING_GROUP, DL1_PIXEL_STATISTICS_GROUP
from ctapipe.io.metadata import read_reference_metadata

from conftest import DARK_FILE, OBSERVATION_RUNS, WAVEFORM_CHUNK_SHIFT, WAVEFORM_CHUNK_SIZE, run_monitoring
from sst1mpipe.io import R0_PIXEL_STATISTICS_GROUP
from sst1mpipe.io.containers import CameraEventType
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.monitoring import PedestalWaveformAggregator, last_chunk_index

TEL_ID = 22
OBSERVATION_FILE, N_PEDESTALS = OBSERVATION_RUNS[TEL_ID]["path"], OBSERVATION_RUNS[TEL_ID]["n_pedestals"]


def read_pedestal_waveforms(path):
    """time, event id and ADC samples (n_events, n_pixels, n_samples) of the pedestal events of a zfits file"""
    with SST1MEventSource(path) as source:
        events = [
            (event.trigger.time, event.index.event_id, np.array(event.r0.tel[TEL_ID].waveform[0], dtype=np.float64))
            for event in source if event.r0.tel[TEL_ID].event_type == CameraEventType.INTERNAL
        ]
    times, event_ids, waveforms = zip(*events, strict=True)
    return Time(times), np.array(event_ids), np.array(waveforms)


@pytest.fixture(scope="module")
def monitoring_file(observation_monitoring):
    """monitoring file of the observation of tel 22, rolling chunks of WAVEFORM_CHUNK_SIZE pedestal events"""
    return observation_monitoring[TEL_ID]


def test_waveform_statistics(monitoring_file):
    """
    statistics of the ADC samples of the pedestal events, in rolling chunks: the 5 last pedestal
    events, every 2 pedestal events
    """
    table = read_table(monitoring_file, f"{R0_PIXEL_STATISTICS_GROUP}/sky_pedestal_waveform/tel_{TEL_ID:03d}")
    times, event_ids, waveforms = read_pedestal_waveforms(OBSERVATION_FILE)
    assert len(times) == N_PEDESTALS == 13

    # overlapping chunks of 5 events shifted by 2 events, the last one ending at the last pedestal
    # event (the last chunk of ctapipe, overlap policy, is the same chunk and is not repeated)
    assert (WAVEFORM_CHUNK_SIZE, WAVEFORM_CHUNK_SHIFT) == (5, 2)
    chunks = [slice(start, start + 5) for start in (0, 2, 4, 6, 8)]
    assert len(table) == len(chunks)
    for row, chunk in zip(table, chunks, strict=True):
        assert row["time_start"] == times[chunk][0] and row["time_end"] == times[chunk][-1]
        assert (row["event_id_start"], row["event_id_end"]) == (event_ids[chunk][0], event_ids[chunk][-1])
        means, stds = waveforms[chunk].mean(axis=-1), waveforms[chunk].std(axis=-1)
        # mean and median of the baseline, mean of the noise of the samples
        np.testing.assert_allclose(row["mean"][0], means.mean(axis=0), rtol=1e-5)
        np.testing.assert_allclose(row["median"][0], np.median(means, axis=0), rtol=1e-5)
        np.testing.assert_allclose(row["std"][0], stds.mean(axis=0), rtol=1e-5)
        assert np.all(row["n_events"] == 5)

    # the dead pixels (low noise) are the outliers, and the chunks are valid
    assert table["std"].shape == (len(chunks), 1, 1296)
    assert np.array_equal(table["outlier_mask"], table["std"] < 0.5 * np.median(table["std"], axis=2, keepdims=True))
    assert np.all(table["is_valid"])
    assert 200 < np.median(table["mean"]) < 600


def test_image_statistics(monitoring_file):
    """statistics of the images of the pedestal events, a single chunk (allow_undersized_tables)"""
    table = read_table(monitoring_file, f"{DL1_PIXEL_STATISTICS_GROUP}/sky_pedestal_image/tel_{TEL_ID:03d}")
    times, _, _ = read_pedestal_waveforms(OBSERVATION_FILE)

    assert len(table) == 1
    assert table["time_start"][0] == times[0] and table["time_end"][0] == times[-1]
    assert table["std"].shape == (1, 1, 1296)
    assert np.all(table["n_events"] == N_PEDESTALS)
    assert not np.any(table["outlier_mask"]) and np.all(table["is_valid"])
    assert 0.1 < np.median(table["std"]) < 3


def test_pointing_table(monitoring_file):
    """pointing of the telescope from its first to its last event"""
    table = read_table(monitoring_file, f"{DL0_TEL_POINTING_GROUP}/tel_{TEL_ID:03d}")
    with SST1MEventSource(OBSERVATION_FILE) as source:
        times = Time([event.trigger.time for event in source])
        pointing, subarray = source.pointing, source.subarray

    assert table["time"][0] == times[0] and table["time"][-1] == times[-1]
    assert np.all(np.diff(table["time"].mjd) > 0)
    # alt/az of the RA/Dec of the run at the location of the telescope
    location = subarray.tel_coords.to_earth_location()[subarray.tel_index_array[TEL_ID]]
    expected = pointing.transform_to(AltAz(obstime=table["time"], location=location))
    filled = SkyCoord(az=table["azimuth"].quantity, alt=table["altitude"].quantity, frame=expected)
    assert np.all(filled.separation(expected) < 1e-6 * u.arcsec)
    assert np.all((74 * u.deg < expected.alt) & (expected.alt < 75 * u.deg))


def test_read_by_ctapipe_monitoring_source(monitoring_file):
    """the monitoring file is read by the ctapipe HDF5MonitoringSource"""
    source = HDF5MonitoringSource(input_files=[monitoring_file])
    assert source.subarray.tel_ids.tolist() == [TEL_ID]
    assert {t.name for t in source.monitoring_types} == {"PIXEL_STATISTICS", "TELESCOPE_POINTINGS"}

    times, _, _ = read_pedestal_waveforms(OBSERVATION_FILE)
    event = ArrayEventContainer()
    event.trigger.time = times[5]
    source.fill_monitoring_container(event)
    monitoring = event.monitoring.tel[TEL_ID]
    image = read_table(monitoring_file, f"{DL1_PIXEL_STATISTICS_GROUP}/sky_pedestal_image/tel_{TEL_ID:03d}")
    np.testing.assert_array_equal(monitoring.camera.pixel_statistics.pedestal_image.std, image["std"][0])
    assert 70 * u.deg < monitoring.pointing.altitude < 80 * u.deg


def test_configuration_and_metadata(monitoring_file):
    assert SubarrayDescription.from_hdf(monitoring_file).tel_ids.tolist() == [TEL_ID]
    scheduling_block = read_table(monitoring_file, "/configuration/observation/scheduling_block")
    assert scheduling_block["sb_type"][0] == SchedulingBlockType.OBSERVATION.value
    observation_block = read_table(monitoring_file, "/configuration/observation/observation_block")
    assert observation_block["target"][0] == "MRK421"
    assert read_reference_metadata(monitoring_file).activity.name == "sst1mpipe-monitoring"


def test_dark_run(tmp_path):
    """dark run: statistics of the dark pedestal events, no pointing"""
    output = tmp_path / "dark.monitoring.h5"
    run_monitoring(DARK_FILE, output)

    with tables.open_file(output) as h5:
        assert DL0_TEL_POINTING_GROUP not in h5
        assert f"{R0_PIXEL_STATISTICS_GROUP}/sky_pedestal_waveform" not in h5
    waveform = read_table(output, f"{R0_PIXEL_STATISTICS_GROUP}/dark_pedestal_waveform/tel_{TEL_ID:03d}")
    image = read_table(output, f"{DL1_PIXEL_STATISTICS_GROUP}/dark_pedestal_image/tel_{TEL_ID:03d}")
    # 225 pedestal events: rolling chunks of the 100 last pedestal events at each pedestal event
    # (126 chunks from the events 0-99 to 125-224, the last chunk ending at the last pedestal
    # event), a single chunk of images
    assert len(waveform) == 126 and len(image) == 1
    assert np.all(waveform["n_events"] == 100)
    with SST1MEventSource(DARK_FILE) as source:
        event_ids = np.array([e.index.event_id for e in source if e.r0.tel[TEL_ID].event_type == CameraEventType.INTERNAL])
    starts = list(range(126))
    np.testing.assert_array_equal(waveform["event_id_start"], event_ids[starts])
    np.testing.assert_array_equal(waveform["event_id_end"], event_ids[np.array(starts) + 99])
    assert np.all(image["n_events"] == 225)
    # no outlier in the dark run (the noise is low in all the pixels)
    assert not np.any(waveform["outlier_mask"]) and np.all(waveform["is_valid"])
    assert np.median(waveform["std"]) < 2.5


def test_chunk_shift_larger_than_chunks(tmp_path):
    """the pedestal events between two chunks would not be used"""
    with pytest.raises(ToolConfigurationError, match="chunk_shift"):
        run_monitoring(
            OBSERVATION_FILE, tmp_path / "monitoring.h5",
            "--ImageStatisticsCalculator.PlainAggregator.SizeChunking.chunk_shift=2000",
        )


def test_last_chunk_index():
    """the statistics at a time are those of the last chunk ended before it"""
    # rolling chunks of 4 pedestal events every 2 pedestal events, pedestal events at 0, 1, ..., 9
    time_start = np.array([0, 2, 4, 6])
    time_end = time_start + 3
    times = [-1, 0, 2.9, 3, 4.5, 5, 6.99, 8, 9, 20]
    # before the end of the first chunk: the first chunk
    assert last_chunk_index(time_start, time_end, times).tolist() == [0, 0, 0, 0, 0, 1, 1, 2, 3, 3]


def test_output_not_overwritten(tmp_path):
    output = tmp_path / "monitoring.h5"
    output.touch()
    with pytest.raises(ToolConfigurationError, match="exists"):
        run_monitoring(OBSERVATION_FILE, output)


def test_pedestal_waveform_aggregator():
    """mean and median of the means, mean of the std, without the masked and invalid values"""
    rng = np.random.default_rng(0)
    data = rng.normal(size=(10, 2, 1, 4))
    data[3, :, 0, 1] = np.nan
    stats = PedestalWaveformAggregator()._compute_stats(data, None)
    valid = np.isfinite(data[:, 0, 0, 1])
    np.testing.assert_allclose(stats.mean[0, 0], data[:, 0, 0, 0].mean())
    np.testing.assert_allclose(stats.median[0, 1], np.median(data[valid, 0, 0, 1]))
    np.testing.assert_allclose(stats.std[0, 2], data[:, 1, 0, 2].mean())
    assert stats.n_events.tolist() == [[10, 9, 10, 10]]

    masked = np.zeros((2, 1, 4), dtype=bool)
    masked[:, 0, 3] = True
    stats = PedestalWaveformAggregator()._compute_stats(data, masked)
    assert np.isnan(stats.mean[0, 3]) and stats.n_events[0, 3] == 0

