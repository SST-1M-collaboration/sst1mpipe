"""
Monitoring data of a SST-1M run, used afterwards by sst1mpipe-process (monitoring source).

The raw data of the run (R0, SST1MEventSource) are read once:

- pointing of the telescopes: altitude and azimuth of the target of the run (TARGET field of
  the file, or SST1MEventSource.pointing_ra/pointing_dec) every pointing_update_interval
  from the first to the last event of the telescope, in the table
  dl0/monitoring/telescope/pointing/tel_XXX (read by the ctapipe PointingInterpolator)
- statistics of the ADC samples of the pedestal events (WaveformStatisticsCalculator), in chunks
  of pedestal events: mean and median of the baseline, std (noise) of the samples, in the table
  r0/monitoring/telescope/calibration/camera/pixel_statistics/<event type>_waveform/tel_XXX
- statistics of the images of the pedestal events (ImageStatisticsCalculator), calibrated by the
  R0R1Calibrator and the CameraCalibrator, in chunks of pedestal events, in the ctapipe pixel
  statistics table
  dl1/monitoring/telescope/calibration/camera/pixel_statistics/<event type>_image/tel_XXX

The event type of the pedestal events is sky_pedestal, dark_pedestal (dark runs) or
electronic_pedestal (bias runs). The statistics are computed as in the
ctapipe-calculate-pixel-statistics tool (aggregation in chunks, outlier detection, validity of
the chunks), and the file is read by the ctapipe HDF5MonitoringSource.
"""
from collections import defaultdict

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz
from astropy.table import QTable, Table
from astropy.time import Time
from ctapipe.calib import CameraCalibrator
from ctapipe.core import Provenance, Tool, ToolConfigurationError
from ctapipe.core.traits import Bool, Float, Path, classes_with_traits, flag
from ctapipe.io import EventSource, HDF5TableWriter, metadata, write_table
from ctapipe.io.datawriter import write_reference_metadata_headers
from ctapipe.io.hdf5dataformat import DL0_TEL_POINTING_GROUP, DL1_PIXEL_STATISTICS_GROUP
from ctapipe.monitoring.aggregator import BaseChunking, StatisticsAggregator
from ctapipe.monitoring.outlier import OutlierDetector
from tqdm import tqdm

from sst1mpipe.calib import R0R1Calibrator
from sst1mpipe.io import R0_PIXEL_STATISTICS_GROUP
from sst1mpipe.io.containers import PEDESTAL_EVENT_TYPES
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.monitoring import (
    ImageStatisticsCalculator,
    WaveformStatisticsCalculator,
    compute_pixel_statistics,
)


class MonitoringTool(Tool):
    """
    Monitoring data of a SST-1M run (raw data): pointing of the telescopes and statistics of
    the ADC samples and of the images of the pedestal events, in chunks of pedestal events
    """

    name = "sst1mpipe-monitoring"
    description = __doc__
    examples = (
        "sst1mpipe-monitoring -i SST1M1_20260120_1179.fits.fz -o SST1M1_20260120_1179.monitoring.h5"
        " --config sst1mpipe/resources/config/sst1mpipe_monitoring_config.json"
    )

    output_path = Path(
        default_value=None,
        allow_none=True,
        directory_ok=False,
        help="Output monitoring file (HDF5)",
    ).tag(config=True)

    overwrite = Bool(default_value=False, help="Overwrite the output file if it exists").tag(config=True)

    progress_bar = Bool(default_value=False, help="Show a progress bar during the event loop").tag(config=True)

    pointing_update_interval = Float(
        default_value=1.0,
        min=1e-3,
        help="Time between two rows of the pointing table of a telescope (s)",
    ).tag(config=True)

    aliases = {
        ("i", "input"): "EventSource.input_url",
        ("o", "output"): "MonitoringTool.output_path",
        ("t", "allowed-tels"): "EventSource.allowed_tels",
        ("m", "max-events"): "EventSource.max_events",
    }

    flags = {
        **flag(
            "overwrite",
            "MonitoringTool.overwrite",
            "overwrite the output file",
            "do not overwrite the output file",
        ),
        **flag(
            "progress",
            "MonitoringTool.progress_bar",
            "show a progress bar during the event loop",
            "do not show a progress bar during the event loop",
        ),
    }

    classes = (
        [R0R1Calibrator, CameraCalibrator, WaveformStatisticsCalculator, ImageStatisticsCalculator]
        + classes_with_traits(EventSource)
        + classes_with_traits(StatisticsAggregator)
        + classes_with_traits(BaseChunking)
        + classes_with_traits(OutlierDetector)
    )

    def setup(self):
        if self.output_path is None:
            raise ToolConfigurationError("The output file is required (-o/--output or MonitoringTool.output_path)")
        if self.output_path.exists() and not self.overwrite:
            raise ToolConfigurationError(f"Output file {self.output_path} exists, use --overwrite to overwrite it")

        self.event_source = self.enter_context(EventSource(parent=self))
        if not isinstance(self.event_source, SST1MEventSource):
            raise ToolConfigurationError(
                f"{self.name} needs the SST-1M raw data (SST1MEventSource), not {type(self.event_source).__name__}"
            )
        subarray = self.event_source.subarray
        self.r0_r1_calibrator = R0R1Calibrator(parent=self, subarray=subarray)
        self.camera_calibrator = CameraCalibrator(parent=self, subarray=subarray)
        self.waveform_calculator = WaveformStatisticsCalculator(parent=self, subarray=subarray)
        self.image_calculator = ImageStatisticsCalculator(parent=self, subarray=subarray)

        # time of the first and last event of each telescope
        self.event_times = {}
        # quantities of the pedestal events of each (telescope, event type), in the event loop
        self.pedestal_waveforms = defaultdict(list)
        self.pedestal_images = defaultdict(list)

    def start(self):
        for event in tqdm(
            self.event_source,
            desc=self.event_source.__class__.__name__,
            total=self.event_source.max_events,
            disable=not self.progress_bar,
        ):
            for tel_id in event.trigger.tels_with_trigger:
                time = event.trigger.tel[tel_id].time
                first, _ = self.event_times.get(tel_id, (time, time))
                self.event_times[tel_id] = (first, time)

            event_type = event.trigger.event_type
            if event_type not in PEDESTAL_EVENT_TYPES:
                continue
            self.add_pedestal_event(event, event_type)

    def add_pedestal_event(self, event, event_type):
        """
        Mean and std of the ADC samples of each pixel (R0), and calibrated image (DL1) of a
        pedestal event. The images are calibrated without the statistics of the ADC samples
        (no voltage drop correction and no dead pixels).
        """
        for tel_id, r0 in event.r0.tel.items():
            samples = np.asarray(r0.waveform, dtype=np.float64)
            moments = np.stack([samples.mean(axis=-1), samples.std(axis=-1)]).astype(np.float32)
            self.pedestal_waveforms[tel_id, event_type].append(
                (event.trigger.tel[tel_id].time, event.index.event_id, moments)
            )
            self.r0_r1_calibrator(event, tel_id)

        self.camera_calibrator(event)
        for tel_id, dl1 in event.dl1.tel.items():
            # shape (n_channels, n_pixels) of the ctapipe pixel statistics
            image = np.asarray(dl1.image, dtype=np.float32)[np.newaxis]
            self.pedestal_images[tel_id, event_type].append((event.trigger.tel[tel_id].time, event.index.event_id, image))

    def finish(self):
        tel_ids = sorted(self.event_times)
        if len(tel_ids) == 0:
            self.log.warning("No event read, no monitoring data written")
            return
        self.write_configuration(tel_ids)
        self.write_pointing_tables()
        self.write_statistics_tables(
            self.pedestal_waveforms, self.waveform_calculator, "waveform", R0_PIXEL_STATISTICS_GROUP,
        )
        self.write_statistics_tables(
            self.pedestal_images, self.image_calculator, "image", DL1_PIXEL_STATISTICS_GROUP,
        )
        Provenance().add_output_file(str(self.output_path), role="Monitoring")
        self.log.info("Monitoring data written in %s", self.output_path)

    def write_configuration(self, tel_ids):
        """
        Description of the telescopes of the run, scheduling and observation blocks and reference
        metadata (as in the files of the ctapipe DataWriter)
        """
        source = self.event_source
        subarray = source.subarray.select_subarray(tel_ids)
        with HDF5TableWriter(self.output_path, mode="w", add_prefix=True, parent=self) as writer:
            subarray.to_hdf(writer.h5file)
            for scheduling_block in source.scheduling_blocks.values():
                writer.write("configuration/observation/scheduling_block", scheduling_block)
            for observation_block in source.observation_blocks.values():
                writer.write("configuration/observation/observation_block", observation_block)
            write_reference_metadata_headers(
                obs_ids=source.obs_ids,
                subarray=subarray,
                writer=writer,
                is_simulation=False,
                data_levels=[],
                contact_info=metadata.Contact(parent=self),
                instrument_info=metadata.Instrument(parent=self),
            )

    def write_pointing_tables(self):
        """
        Altitude and azimuth of the pointing of each telescope every pointing_update_interval,
        from its first to its last event
        """
        pointing = self.event_source.pointing
        if pointing is None:
            self.log.warning("No pointing for the run (target %r): no pointing table", self.event_source.target)
            return
        subarray = self.event_source.subarray
        locations = subarray.tel_coords.to_earth_location()
        for tel_id, (first, last) in self.event_times.items():
            duration = (last - first).to_value(u.s)
            offsets = np.unique(np.append(np.arange(0, duration, self.pointing_update_interval), duration))
            times = first + offsets * u.s
            location = locations[subarray.tel_index_array[tel_id]]
            altaz = pointing.transform_to(AltAz(obstime=times, location=location))
            table = QTable(dict(time=times, azimuth=altaz.az.to(u.rad), altitude=altaz.alt.to(u.rad)))
            write_table(table, self.output_path, f"{DL0_TEL_POINTING_GROUP}/tel_{tel_id:03d}")
            self.log.info("Pointing of telescope %d: %d rows", tel_id, len(table))

    def write_statistics_tables(self, pedestal_events, calculator, column, group):
        """
        Statistics of the quantities ``column`` of the pedestal events of each telescope and
        event type, in the tables ``group/<event type>_<column>/tel_XXX``
        """
        for (tel_id, event_type), rows in sorted(pedestal_events.items(), key=lambda item: (item[0][0], item[0][1].value)):
            times, event_ids, values = zip(*rows, strict=True)
            table = Table(dict(time=Time(times), event_id=np.array(event_ids), **{column: np.array(values)}))
            name = f"{event_type.name.lower()}_{column}"
            try:
                statistics = compute_pixel_statistics(calculator, table, tel_id, column)
            except ValueError as error:
                self.log.warning("No statistics %s of telescope %d (%d events): %s", name, tel_id, len(table), error)
                continue
            write_table(statistics, self.output_path, f"{group}/{name}/tel_{tel_id:03d}")
            self.log.info(
                "Statistics %s of telescope %d: %d pedestal events, %d chunks (%d valid)",
                name, tel_id, len(table), len(statistics), np.count_nonzero(statistics["is_valid"]),
            )


def main():
    tool = MonitoringTool()
    tool.run()


if __name__ == "__main__":
    main()
