"""
Statistics of the pedestal events of the SST-1M telescopes, in chunks of events.

The statistics are computed with the ctapipe `~ctapipe.monitoring.PixelStatisticsCalculator`
(aggregation in chunks of events, outlier detection, validity of the chunks), as in the
``ctapipe-calculate-pixel-statistics`` tool:

- `WaveformStatisticsCalculator`: statistics of the ADC samples of the pedestal events (R0),
  aggregated by the `PedestalWaveformAggregator`
- `ImageStatisticsCalculator`: statistics of the calibrated images of the pedestal events (DL1)

The two calculators are distinct classes, so that they are configured separately,
e.g. ``{"WaveformStatisticsCalculator": {"PedestalWaveformAggregator": {"SizeChunking": {"chunk_size": 100}}}}``.
"""
import numpy as np
from astropy.table import vstack
from ctapipe.containers import ChunkStatisticsContainer
from ctapipe.core.traits import ComponentName, TelescopeParameter
from ctapipe.monitoring.aggregator import StatisticsAggregator
from ctapipe.monitoring.calculator import PixelStatisticsCalculator

__all__ = [
    "PedestalWaveformAggregator",
    "WaveformStatisticsCalculator",
    "ImageStatisticsCalculator",
    "compute_pixel_statistics",
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
    Statistics of the column ``col_name`` of ``table`` in chunks of events, as in the
    ``ctapipe-calculate-pixel-statistics`` tool: a first pass with non-overlapping chunks, and
    a second pass with overlapping chunks (``chunk_shift``) over the faulty chunks.

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
        statistics of the chunks: time_start, time_end, event_id_start, event_id_end,
        n_events, mean, median, std, outlier_mask, is_valid
    """
    statistics = calculator.first_pass(table=table, tel_id=tel_id, col_name=col_name)
    aggregator = calculator.stats_aggregators[calculator.stats_aggregator_type.tel[tel_id]]
    chunk_shift = getattr(aggregator.chunking, "chunk_shift", None)
    has_shift = chunk_shift is not None and np.any(chunk_shift > 0)
    if has_shift and np.any(~statistics["is_valid"].data):
        second_pass = calculator.second_pass(
            table=table, valid_chunks=statistics["is_valid"].data, tel_id=tel_id, col_name=col_name,
        )
        statistics = vstack([statistics, second_pass])
        statistics.sort(["time_start"])
    return statistics
