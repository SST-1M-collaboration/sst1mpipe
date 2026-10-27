"""
Monitoring source of the monitoring files of sst1mpipe-monitoring, used by sst1mpipe-process
(``--monitoring-source SST1MMonitoringSource --SST1MMonitoringSource.input_files=<file>``).
"""
import astropy.units as u
from astropy.time import Time
import numpy as np
import tables
from ctapipe.containers import StatisticsContainer, TelescopePointingContainer
from ctapipe.io import HDF5MonitoringSource, read_table
from ctapipe.io.monitoringtypes import MonitoringType
from ctapipe.monitoring import StatisticsInterpolator

from sst1mpipe.io.io import R0_PIXEL_STATISTICS_GROUP
from sst1mpipe.monitoring import last_chunk_index

__all__ = ["LastChunkStatisticsInterpolator", "SST1MMonitoringSource", "statistics_container", "waveform_statistics_interpolator"]


class LastChunkStatisticsInterpolator(StatisticsInterpolator):
    """
    Statistics of the last chunk of pedestal events ended at or before the requested time
    (`~sst1mpipe.monitoring.last_chunk_index`): with rolling chunks, the statistics of the last
    ``chunk_size`` pedestal events before the event. The ctapipe `~ctapipe.monitoring.ChunkInterpolator`
    takes instead the last chunk started before the time, which contains the pedestal events
    after it. The times before the end of the first chunk use the first chunk. Out of the time
    range of the chunks (within ``timestamp_tolerance``), the statistics are NaN.
    """

    def add_table(self, tel_id, input_table):
        """As the ctapipe ChunkInterpolator, but the start and end of the chunks are kept as astropy times"""
        self._check_tables(input_table)
        table = input_table.copy()
        table.sort("time_start")
        self.time_start[tel_id] = Time(table["time_start"])
        self.time_end[tel_id] = Time(table["time_end"])
        self.values[tel_id] = {column: table[column] for column in self.columns}

    def _interpolate_chunk(self, tel_id, column, time, timestamp_tolerance=0 * u.s):
        time_start, time_end = self.time_start[tel_id], self.time_end[tel_id]
        values = self.values[tel_id][column]
        times = time.reshape(-1)
        index = last_chunk_index(time_end, times)
        # astropy times: compared in a common time scale
        in_range = (times >= time_start[0] - timestamp_tolerance) & (times <= time_end.max() + timestamp_tolerance)
        no_value = np.nan if np.isscalar(values[0]) else np.full_like(values[0], np.nan)
        interpolated = [values[i] if ok else no_value for i, ok in zip(index, in_range, strict=True)]
        return interpolated[0] if len(interpolated) == 1 else interpolated


class SST1MMonitoringSource(HDF5MonitoringSource):
    """
    Monitoring data of the SST-1M runs (files of sst1mpipe-monitoring), interpolated at the
    time of the events in their monitoring containers (``event.monitoring.tel[tel_id]``),
    as by the ctapipe `~ctapipe.io.HDF5MonitoringSource`, with:

    - the statistics of the ADC samples of the pedestal events (R0 tables ``<event type>_waveform``)
      in ``camera.pixel_statistics.pedestal_waveform``, used by the R0R1Calibrator
      (voltage drop correction and dead pixels)
    - the statistics of the images of the pedestal events of channel 0 in
      ``camera.pixel_statistics.pedestal_image``, of shape (n_pixels,) as the images (DigiCam
      has a single channel), for the NSBImageCleaner

    The statistics at the time of the event are those of the last chunk of pedestal events ended
    before it (`LastChunkStatisticsInterpolator`): with the rolling chunks of sst1mpipe-monitoring,
    the statistics of the last pedestal events. Out of the time range of the chunks, within
    ``timestamp_tolerance``, the statistics are not filled (``std`` is None). The outlier pixels of
    the statistics (e.g. the dead pixels) are NaN.
    Out of the time range of the pointing table, the pointing is NaN.
    """

    def _process_single_file(self, file):
        super()._process_single_file(file)
        self._process_image_statistics()
        self._process_waveform_statistics(file)

    def _process_image_statistics(self):
        """Statistics of the images of the pedestal events: last chunk ended before the event"""
        for name in list(getattr(self, "pixel_stats_dict", {})):
            interpolator = LastChunkStatisticsInterpolator()
            for tel_id in self.subarray.tel_ids:
                # the outliers are already NaN (HDF5MonitoringSource)
                interpolator.add_table(tel_id, self._pixel_statistics[tel_id][name])
            self.pixel_stats_dict[name] = interpolator

    def _process_waveform_statistics(self, file):
        """Statistics of the ADC samples of the pedestal events (R0 pixel statistics tables)"""
        if not hasattr(self, "_waveform_interpolators"):
            self._waveform_interpolators = {}
        with tables.open_file(file) as h5file:
            if R0_PIXEL_STATISTICS_GROUP not in h5file:
                return
            names = [group._v_name for group in h5file.get_node(R0_PIXEL_STATISTICS_GROUP)._f_iter_nodes("Group")]
        for name in names:
            self._waveform_interpolators[name] = waveform_statistics_interpolator({
                tel_id: read_table(file, f"{R0_PIXEL_STATISTICS_GROUP}/{name}/tel_{tel_id:03d}")
                for tel_id in self.subarray.tel_ids
            })

    def fill_monitoring_container(self, event):
        """
        Fill the monitoring containers of the telescopes of the monitoring files, at the
        time of the event
        """
        time = event.trigger.time
        for tel_id in self.subarray.tel_ids:
            monitoring = event.monitoring.tel[tel_id]
            statistics = monitoring.camera.pixel_statistics
            if self.has_pixel_statistics:
                for name in self.pixel_stats_dict:
                    values = self.get_values(MonitoringType.PIXEL_STATISTICS, time=time, tel_id=tel_id, subtype=name)
                    field = "pedestal_image" if "pedestal_image" in name else name
                    # statistics of channel 0: the images are gain selected
                    statistics[field] = statistics_container(values, channel=0)
                for interpolator in getattr(self, "_waveform_interpolators", {}).values():
                    if "pedestal_waveform" in statistics.fields:
                        values = interpolator(tel_id, time, self.timestamp_tolerance)
                        statistics["pedestal_waveform"] = statistics_container(values)
            if self.has_camera_coefficients:
                monitoring.camera.coefficients = self.get_camera_monitoring_container(tel_id, time).coefficients
            if self.has_pointings:
                try:
                    monitoring.pointing = self.get_telescope_pointing_container(tel_id, time)
                except ValueError:
                    # out of the time range of the pointing table
                    monitoring.pointing = TelescopePointingContainer()


def waveform_statistics_interpolator(tables):
    """
    Interpolator of the statistics of the ADC samples of the pedestal events (R0 tables
    ``<event type>_waveform`` of sst1mpipe-monitoring) of the telescopes, ``tables`` by telescope id:
    the last chunk ended before the time (`LastChunkStatisticsInterpolator`), the outliers (e.g. the
    dead pixels) NaN, as in the HDF5MonitoringSource. Used by SST1MMonitoringSource and by
    sst1mpipe-monitoring to calibrate the images of the pedestal events.
    """
    interpolator = LastChunkStatisticsInterpolator()
    for tel_id, table in tables.items():
        table = table.copy()
        for column in ("mean", "median", "std"):
            table[column][table["outlier_mask"].data] = np.nan
        interpolator.add_table(tel_id, table)
    return interpolator


def statistics_container(values, channel=None):
    """
    StatisticsContainer of interpolated statistics (mean, median, std), of a channel if given.
    Not filled (None) if there are no statistics at this time (all NaN).
    """
    if np.all(np.isnan(values["std"])):
        return StatisticsContainer()
    if channel is not None:
        values = {column: values[column][channel] for column in ("mean", "median", "std")}
    return StatisticsContainer(mean=values["mean"], median=values["median"], std=values["std"])
