"""
Monitoring source of the monitoring files of sst1mpipe-monitoring, used by sst1mpipe-process
(``--monitoring-source SST1MMonitoringSource --SST1MMonitoringSource.input_files=<file>``).
"""
import numpy as np
import tables
from ctapipe.containers import StatisticsContainer, TelescopePointingContainer
from ctapipe.io import HDF5MonitoringSource, read_table
from ctapipe.io.monitoringtypes import MonitoringType
from ctapipe.monitoring import StatisticsInterpolator

from sst1mpipe.io.io import R0_PIXEL_STATISTICS_GROUP

__all__ = ["SST1MMonitoringSource"]


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

    The statistics are interpolated with the chunks of pedestal events containing the time of the
    event, within ``timestamp_tolerance``. Without such a chunk, the statistics are not filled
    (``std`` is None). The outlier pixels of the statistics (e.g. the dead pixels) are NaN.
    Out of the time range of the pointing table, the pointing is NaN.
    """

    def _process_single_file(self, file):
        super()._process_single_file(file)
        self._process_waveform_statistics(file)

    def _process_waveform_statistics(self, file):
        """Statistics of the ADC samples of the pedestal events (R0 pixel statistics tables)"""
        if not hasattr(self, "_waveform_interpolators"):
            self._waveform_interpolators = {}
        with tables.open_file(file) as h5file:
            if R0_PIXEL_STATISTICS_GROUP not in h5file:
                return
            names = [group._v_name for group in h5file.get_node(R0_PIXEL_STATISTICS_GROUP)._f_iter_nodes("Group")]
        for name in names:
            interpolator = StatisticsInterpolator()
            for tel_id in self.subarray.tel_ids:
                table = read_table(file, f"{R0_PIXEL_STATISTICS_GROUP}/{name}/tel_{tel_id:03d}")
                # the outliers (e.g. the dead pixels) are NaN, as in the HDF5MonitoringSource
                for column in ("mean", "median", "std"):
                    table[column][table["outlier_mask"].data] = np.nan
                interpolator.add_table(tel_id, table)
            self._waveform_interpolators[name] = interpolator

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
