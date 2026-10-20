"""
Statistics of the pedestal events of the SST-1M telescopes, in chunks of events.

The statistics are computed with the ctapipe `~ctapipe.monitoring.PixelStatisticsCalculator`
(aggregation in chunks of events, outlier detection, validity of the chunks), as in the
``ctapipe-calculate-pixel-statistics`` tool:

- `WaveformStatisticsCalculator`: statistics of the ADC samples of the pedestal events (R0),
  aggregated by the `PedestalWaveformAggregator`
- `ImageStatisticsCalculator`: statistics of the calibrated images of the pedestal events (DL1)

The chunks are rolling with ``chunk_shift``: a new chunk every ``chunk_shift`` pedestal events, each
of the last ``chunk_size`` pedestal events (`compute_pixel_statistics`). At the time of an event, the
statistics are those of the last chunk ended before it (`last_chunk_index`).

The two calculators are distinct classes, so that they are configured separately,
e.g. ``{"WaveformStatisticsCalculator": {"PedestalWaveformAggregator": {"SizeChunking": {"chunk_size": 100}}}}``.
"""
import numpy as np
from ctapipe.containers import ChunkStatisticsContainer
from ctapipe.core.traits import ComponentName, TelescopeParameter
from ctapipe.monitoring.aggregator import SizeChunking, StatisticsAggregator, TimeChunking
from ctapipe.monitoring.calculator import PixelStatisticsCalculator

__all__ = [
    "PedestalWaveformAggregator",
    "WaveformStatisticsCalculator",
    "ImageStatisticsCalculator",
    "check_chunk_shift",
    "compute_pixel_statistics",
    "last_chunk_index",
]


class PedestalWaveformAggregator(StatisticsAggregator):
    """
    Statistics of the ADC samples of the pedestal events, from the mean and the std of the
    samples of each pixel of each event: ``data`` of shape (n_events, 2, n_channels, n_pixels),
    ``data[:, 0]`` the mean and ``data[:, 1]`` the std of the samples.

    - mean, median: mean and median over the events of the mean of the samples (baseline)
    - std: mean over the events of the std of the samples (noise, used for the voltage drop
      correction and the dead pixels)
    """

    def _compute_stats(self, data, masked_elements_of_sample) -> ChunkStatisticsContainer:
        mask = None
        if masked_elements_of_sample is not None:
            mask = np.broadcast_to(masked_elements_of_sample, np.shape(data))
        data = np.ma.masked_invalid(np.ma.array(data, mask=mask))
        means, stds = data[:, 0], data[:, 1]
        return ChunkStatisticsContainer(
            n_events=np.count_nonzero(~np.ma.getmaskarray(means), axis=0),
            mean=np.ma.filled(np.ma.mean(means, axis=0), np.nan),
            median=np.ma.filled(np.ma.median(means, axis=0), np.nan),
            std=np.ma.filled(np.ma.mean(stds, axis=0), np.nan),
        )


class WaveformStatisticsCalculator(PixelStatisticsCalculator):
    """
    Statistics of the ADC samples of the pedestal events (R0), see `PedestalWaveformAggregator`
    """

    stats_aggregator_type = TelescopeParameter(
        trait=ComponentName(StatisticsAggregator, default_value="PedestalWaveformAggregator"),
        default_value="PedestalWaveformAggregator",
        help="Name of the StatisticsAggregator subclass to be used.",
    ).tag(config=True)


class ImageStatisticsCalculator(PixelStatisticsCalculator):
    """Statistics of the calibrated images (p.e.) of the pedestal events (DL1)"""


def compute_pixel_statistics(calculator, table, tel_id, col_name):
    """
    Statistics of the column ``col_name`` of ``table`` in chunks of events, with the outlier
    detection and the validity of the chunks of the ``ctapipe-calculate-pixel-statistics`` tool.

    The chunks follow the chunking of the aggregator: with ``chunk_shift``, the chunks are
    rolling, i.e. a new chunk every ``chunk_shift`` events (or duration), each of ``chunk_size``
    events (or ``chunk_duration``), overlapping the previous ones; without ``chunk_shift``, the
    chunks do not overlap. The last chunk ends at the last event (``last_chunk_policy``
    ``overlap``). Unlike the first pass of ctapipe, ``chunk_shift`` is used for all the
    chunks.

    Parameters
    ----------
    calculator: ctapipe.monitoring.PixelStatisticsCalculator
    table: astropy.table.Table
        events of the telescope, with the columns time, event_id and ``col_name``
        (shape (n_events, ...))
    tel_id: int
    col_name: str

    Returns
    -------
    astropy.table.Table
        statistics of the chunks, sorted by time: time_start, time_end, event_id_start,
        event_id_end, n_events, mean, median, std, outlier_mask, is_valid
    """
    aggregator = calculator.stats_aggregators[calculator.stats_aggregator_type.tel[tel_id]]
    statistics = aggregator(table=table, col_name=col_name)
    # with a shift smaller than the chunks, the last chunk of ctapipe (overlap policy) can be
    # the previous chunk again
    duplicate = np.zeros(len(statistics), dtype=bool)
    duplicate[1:] = (statistics["event_id_start"][1:] == statistics["event_id_start"][:-1]) & (
        statistics["event_id_end"][1:] == statistics["event_id_end"][:-1]
    )
    statistics = statistics[~duplicate]
    calculator._find_and_append_outliers(statistics)
    statistics["is_valid"] = calculator._get_valid_chunks(statistics["outlier_mask"])
    return statistics


def check_chunk_shift(calculator):
    """
    Raise a ValueError if the chunks of the aggregators of ``calculator`` are shifted by more than
    their size: the pedestal events between two chunks would not be used
    """
    for name, aggregator in calculator.stats_aggregators.items():
        chunking = aggregator.chunking
        if isinstance(chunking, SizeChunking):
            size, shift = chunking.chunk_size, chunking.chunk_shift
        elif isinstance(chunking, TimeChunking):
            size, shift = chunking.chunk_duration, chunking.chunk_shift
        else:
            continue
        if isinstance(chunking, SizeChunking) and shift is not None and shift < 1:
            raise ValueError(f"{type(calculator).__name__}.{name}.SizeChunking: chunk_shift ({shift}) must be at least 1")
        if size is not None and shift is not None and shift > size:
            raise ValueError(
                f"{type(calculator).__name__}.{name}.{type(chunking).__name__}: chunk_shift ({shift}) "
                f"larger than the chunks ({size})"
            )


def last_chunk_index(time_start, time_end, times):
    """
    Index of the chunk of pedestal events to use at each time: the last chunk ended at or before
    the time, i.e. the statistics of the last pedestal events before it. The times before the
    end of the first chunk use the first chunk.

    Parameters
    ----------
    time_start, time_end: array
        start and end of the chunks (MJD), sorted by time
    times: array
        times (MJD)

    Returns
    -------
    numpy.ndarray of int
        index of the chunk at each time
    """
    # the end of the chunks increases with their start, except if the chunks have different sizes
    time_end = np.maximum.accumulate(np.asarray(time_end))
    index = np.searchsorted(time_end, np.asarray(times), side="right") - 1
    return np.clip(index, 0, len(time_start) - 1)
