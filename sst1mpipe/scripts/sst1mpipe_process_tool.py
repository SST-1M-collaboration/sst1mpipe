import os
from collections import Counter

import astropy.units as u
import numpy as np
from astropy.table import QTable
from astropy.time import Time
from tqdm import tqdm
from ctapipe.calib import CameraCalibrator
from ctapipe.containers import EventType, SchedulingBlockType, TelEventIndexContainer
from ctapipe.core import Tool
from ctapipe.core.traits import Bool, Integer, List, UseEnum, flag
from ctapipe.image import ImageProcessor
from ctapipe.io import EventSource, DataWriter, SimTelEventSource, write_table
from ctapipe.reco import ShowerProcessor

from sst1mpipe.calib import R0R1Calibrator, ImageSaturationCorrector
from sst1mpipe.io import write_dl1_info
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils.monitoring_pedestals import DL1PedestalMonitor, R0PedestalMonitor
from sst1mpipe.utils.cleaning import DBSCANImageCleaner, TimeDBSCANImageCleaner
from sst1mpipe.io.zmq_event_source import ZMQEventSource


class ProcessorTool(Tool):
    """
    Process data from lower-data levels up to DL1 including image
    extraction and optionally image parameterization.
    This implementation is based on the ctapipe.tool.ProcessorTool

    For the SST-1M raw data (R0, SST1MEventSource), the R0 -> R1 calibration is done
    by the R0R1Calibrator. The statistics of the ADC samples of the pedestal events,
    used for the voltage drop correction and the dead pixels, are computed in a sliding
    window by the R0PedestalMonitor: these corrections are applied once the first
    pedestal event of the telescope is read.
    In the same way, the statistics of the calibrated images of the pedestal events are
    computed by the DL1PedestalMonitor, for the NSBImageCleaner which raises the picture
    threshold of the pixels with a high pedestal std.
    The charges and peak times of the saturated pixels are corrected by the ImageSaturationCorrector.
    For the simulations (SimTelEventSource), the R1 waveforms are corrected for the
    PDE drop by the R0R1Calibrator (pde_drop_factor).
    Only the events of the runs with an allowed scheduling block type (allowed_sb_types)
    are processed, e.g. the observations but not the transitions between two wobbles.

    The pedestal events are used for the monitoring only: they are not written as events.
    The statistics of their ADC samples and images are written every
    pedestal_monitoring_interval pedestal events in the tables
    r0/monitoring/telescope/pedestal and dl1/monitoring/telescope/pedestal.
    The pointing of the telescopes is written in dl0/monitoring/telescope/pointing/tel_XXX.
    The shower geometry is reconstructed by the ShowerProcessor if the DL2 is written
    (DataWriter.write_dl2).
    """

    name = 'sst1mpipe-process'
    description = __doc__
    examples = ("sst1mpipe-process -i mysim.simtel.gz -o events.dl1.h5",
                "sst1mpipe-process -i tcp://localhost:24593 -o events.dl1.h5 "
                "--config sst1mpipe/resources/config/sst1mpipe_rta_config.json --log-level INFO")

    progress_bar = Bool(
        help="show progress bar during processing", default_value=False
    ).tag(config=True)

    wobble_in_output_name = Bool(
        help=(
            "Add the wobble of the observation (e.g. W1, read from the input file)"
            " to the name of the output file, e.g. events_W1.dl1.h5"
        ),
        default_value=True,
    ).tag(config=True)

    allowed_sb_types = List(
        UseEnum(SchedulingBlockType),
        default_value=[SchedulingBlockType.OBSERVATION],
        help=(
            "Scheduling block types of the runs whose events are processed, e.g. OBSERVATION,"
            " CALIBRATION (dark runs), ENGINEERING, UNKNOWN (transitions between two wobbles)."
            " The events of the other runs are skipped. The events without scheduling block"
            " (e.g. of a ZMQ stream) are processed."
        ),
    ).tag(config=True)

    pedestal_monitoring_interval = Integer(
        default_value=20,
        min=1,
        help=(
            "Number of pedestal events of a telescope between two rows of the pedestal"
            " monitoring tables"
        ),
    ).tag(config=True)

    aliases = {
        ("i", "input"): "EventSource.input_url",
        ("o", "output"): "DataWriter.output_path",
        ("t", "allowed-tels"): "EventSource.allowed_tels",
        ("m", "max-events"): "EventSource.max_events",
    }

    flags = {

        **flag(
            "progress",
            "ProcessorTool.progress_bar",
            "show a progress bar during event processing",
            "don't show a progress bar during event processing",)
    }

    classes = [
        DBSCANImageCleaner, TimeDBSCANImageCleaner, R0R1Calibrator, R0PedestalMonitor, DL1PedestalMonitor,
        ImageSaturationCorrector, ShowerProcessor,
    ]

    def setup(self):

        if ZMQEventSource.is_compatible(self.config.EventSource.input_url): # temporary fix since tcp:// url is not accepted by EventSource
            self.event_source = self.enter_context(ZMQEventSource(parent=self))
        else:
            self.event_source = self.enter_context(EventSource(parent=self))
        # R0 -> R1 calibration of the SST-1M raw data, PDE drop correction of the simulations.
        # The pedestal statistics are computed for the SST-1M data, from its pedestal events
        self.r0_pedestal_monitor = None
        self.dl1_pedestal_monitor = None
        self.r0_r1_calibrator = None
        self.image_saturation_corrector = None
        subarray = self.event_source.subarray
        # the SST-1M raw data (R0) of the files or of the ZMQ stream (DigiCam camera events);
        # for the R1 events of a stream, the R0 -> R1 calibration does nothing
        if isinstance(self.event_source, SST1MEventSource | ZMQEventSource):
            self.r0_pedestal_monitor = R0PedestalMonitor(parent=self, subarray=subarray)
            self.image_saturation_corrector = ImageSaturationCorrector(parent=self, subarray=subarray)
            self.dl1_pedestal_monitor = DL1PedestalMonitor(parent=self, subarray=subarray)
            self.r0_r1_calibrator = R0R1Calibrator(parent=self, subarray=subarray)
        elif isinstance(self.event_source, SimTelEventSource):
            self.r0_r1_calibrator = R0R1Calibrator(parent=self, subarray=subarray)
        self.camera_calibrator = CameraCalibrator(parent=self, subarray=self.event_source.subarray)
        self.image_processor = ImageProcessor(parent=self, subarray=self.event_source.subarray)
        # the writer is closed in finish(), to read back the output file. If the processing
        # fails before, it is closed when the tool exits.
        self.writer = DataWriter(event_source=self.event_source, parent=self)
        # reconstruction of the shower geometry, written in the DL2
        self.shower_processor = None
        if self.writer.write_dl2:
            self.shower_processor = ShowerProcessor(parent=self, subarray=self.event_source.subarray)
        self._sb_types = {}
        self.n_skipped_events = 0
        # processing summary, counted in the event loop
        self.n_triggered = Counter()
        self.n_pedestal = 0
        self.n_survived_pedestals = 0
        # pointing of each telescope (time, azimuth, altitude): the rows of the pointing table
        # and the pointing of the last event
        self._pointing_rows = {}
        self._last_pointing = {}
        self._writer_closed = False
        self._exit_stack.callback(self._close_writer)

    def _close_writer(self):
        if not self._writer_closed:
            self._writer_closed = True
            self.writer.finish()

    def start(self):

        for event in tqdm(
            self.event_source,
            desc=self.event_source.__class__.__name__,
            total=self.event_source.max_events,
            disable=not self.progress_bar,
        ):
            if not self.is_allowed(event):
                self.n_skipped_events += 1
                continue
            if self.r0_r1_calibrator is not None:
                self.calibrate_r0_r1(event)
            self.camera_calibrator(event)
            # the saturated pixels are corrected with the R0 waveforms
            if self.image_saturation_corrector is not None and len(event.r0.tel) > 0:
                self.image_saturation_corrector(event)
            if self.dl1_pedestal_monitor is not None:
                self.fill_dl1_pedestal_monitoring(event)
            self.image_processor(event)
            if self.dl1_pedestal_monitor is not None:
                self.add_dl1_pedestal(event)
            self.add_pointing(event)
            self.n_triggered.update(event.trigger.tels_with_trigger)

            # the pedestal events are used for the monitoring only
            if event.trigger.event_type == EventType.SKY_PEDESTAL:
                self.n_pedestal += 1
                self.n_survived_pedestals += any(
                    np.isfinite(dl1.parameters.hillas.intensity) for dl1 in event.dl1.tel.values()
                )
                self.write_pedestal_monitoring(event)
                continue

            if self.shower_processor is not None:
                self.shower_processor(event)
            self.writer(event)

    def scheduling_block_type(self, obs_id):
        """Type of the scheduling block of the observation block, None if unknown"""
        if obs_id not in self._sb_types:
            observation_block = self.event_source.observation_blocks.get(obs_id)
            scheduling_block = None
            if observation_block is not None:
                scheduling_block = self.event_source.scheduling_blocks.get(int(observation_block.sb_id))
            self._sb_types[obs_id] = None if scheduling_block is None else scheduling_block.sb_type
            if self._sb_types[obs_id] is not None and self._sb_types[obs_id] not in self.allowed_sb_types:
                self.log.warning(
                    "Events of obs_id %d skipped: scheduling block of type %s, allowed types: %s",
                    obs_id, self._sb_types[obs_id].name, [t.name for t in self.allowed_sb_types],
                )
        return self._sb_types[obs_id]

    def is_allowed(self, event):
        """True if the event belongs to a run with an allowed scheduling block type (or without one)"""
        sb_type = self.scheduling_block_type(event.index.obs_id)
        return sb_type is None or sb_type in self.allowed_sb_types

    def calibrate_r0_r1(self, event):
        """R0 -> R1 calibration, with the pedestal statistics of the sliding window"""
        if self.r0_pedestal_monitor is None:
            # simulation: PDE drop correction of the R1 waveforms
            self.r0_r1_calibrator(event)
            return
        for tel_id in event.r0.tel:
            if event.trigger.event_type == EventType.SKY_PEDESTAL:
                self.r0_pedestal_monitor.add_event(event, tel_id)
            self.r0_pedestal_monitor.fill_monitoring(event, tel_id)
            self.r0_r1_calibrator(event, tel_id)

    def fill_dl1_pedestal_monitoring(self, event):
        """
        Statistics of the images of the pedestal events in event.mon.tel[tel_id].pedestal,
        used by the NSBImageCleaner (nothing is filled before the first pedestal event)
        """
        for tel_id in event.dl1.tel:
            self.dl1_pedestal_monitor.fill_monitoring(event, tel_id)

    def add_dl1_pedestal(self, event):
        """Add the image of a pedestal event to the sliding window, after its cleaning"""
        if event.trigger.event_type == EventType.SKY_PEDESTAL:
            for tel_id in event.dl1.tel:
                self.dl1_pedestal_monitor.add_event(event, tel_id)

    def add_pointing(self, event):
        """
        Pointing of the telescopes (event.pointing.tel) for the tables
        dl0/monitoring/telescope/pointing/tel_XXX read by the ctapipe PointingInterpolator:
        a row each time the alt/az of the pointing is computed by the event source,
        and a last row at the time of the last event (see write_pointing_tables)
        """
        # the pointing of the simulations is written by the DataWriter
        if self.event_source.is_simulation:
            return
        for tel_id in event.trigger.tels_with_trigger:
            pointing = event.pointing.tel[tel_id]
            if not np.isfinite(pointing.altitude):
                continue
            row = (event.trigger.time, pointing.azimuth, pointing.altitude)
            rows = self._pointing_rows.setdefault(tel_id, [])
            if not rows or row[1:] != rows[-1][1:]:
                rows.append(row)
            self._last_pointing[tel_id] = row

    def write_pointing_tables(self, output_path):
        """Pointing tables of the telescopes, written after the events"""
        for tel_id, rows in self._pointing_rows.items():
            if self._last_pointing[tel_id] is not rows[-1]:
                rows.append(self._last_pointing[tel_id])
            time, azimuth, altitude = zip(*rows, strict=True)
            table = QTable(dict(time=Time(time), azimuth=u.Quantity(azimuth), altitude=u.Quantity(altitude)))
            write_table(table, output_path, f"/dl0/monitoring/telescope/pointing/tel_{tel_id:03d}")

    def write_pedestal_monitoring(self, event):
        """
        Statistics of the pedestal events (ADC samples and images) in the pedestal monitoring
        tables, every pedestal_monitoring_interval pedestal events of the telescope
        """
        if self.r0_pedestal_monitor is None:
            return
        for tel_id in event.r0.tel:
            if self.r0_pedestal_monitor.processed_events[tel_id] % self.pedestal_monitoring_interval != 0:
                continue
            index = TelEventIndexContainer(obs_id=event.index.obs_id, event_id=event.index.event_id, tel_id=tel_id)
            monitoring = event.mon.tel[tel_id]
            self.writer._writer.write("r0/monitoring/telescope/pedestal", [index, monitoring.r0])
            # statistics of the images including the image of this event
            self.dl1_pedestal_monitor.fill_monitoring(event, tel_id)
            if monitoring.pedestal.charge_std is not None:
                self.writer._writer.write("dl1/monitoring/telescope/pedestal", [index, monitoring.pedestal])

    def finish(self):

        if self.event_source.is_simulation:
            self.writer.write_simulated_shower_distributions(self.event_source.simulated_shower_distributions)
        self._close_writer()
        output_path = self.writer.output_path
        self.write_pointing_tables(output_path)
        if self.n_skipped_events > 0:
            self.log.warning("%d events skipped (scheduling block type not allowed)", self.n_skipped_events)

        n_triggered = {int(tel_id): self.n_triggered[tel_id] for tel_id in self.event_source.subarray.tel_ids}
        for tel_id, n in n_triggered.items():
            self.log.info("Number of triggered events of telescope %d: %d", tel_id, n)
        self.log.info("Number of pedestal events: %d", self.n_pedestal)
        if self.n_pedestal > 0:
            self.log.info(
                "Fraction of pedestal events that survived cleaning: %f",
                self.n_survived_pedestals / self.n_pedestal,
            )

        # observation information, from the event source (SST1MEventSource)
        source = self.event_source
        target = getattr(source, "target", None)
        wobble = getattr(source, "wobble", None)
        pointing = getattr(source, "pointing", None)
        pointing_manual = getattr(source, "pointing_manual", False)
        self.log.info("Target: %s, wobble: %s, pointing: %s (manual: %s)", target, wobble, pointing, pointing_manual)

        calibration_files, window_files = None, None
        if self.r0_pedestal_monitor is not None:  # observed data
            tel_ids = [tel_id for tel_id, n in n_triggered.items() if n > 0]
            calibration_files = ",".join(str(self.r0_r1_calibrator.calibration_file_path(t)) for t in tel_ids)
            window_files = ",".join(str(self.r0_r1_calibrator.window_transmittance_file_path(t)) for t in tel_ids)

        n_triggered = list(n_triggered.values()) + [0, 0]
        write_dl1_info(output_path, dict(
            calib_file=calibration_files,
            window_file=window_files,
            target=target,
            wobble=wobble,
            ra=None if pointing is None else pointing.ra.deg,
            dec=None if pointing is None else pointing.dec.deg,
            manual_coords=pointing_manual,
            n_saturated=(
                None if self.image_saturation_corrector is None
                else sum(self.image_saturation_corrector.n_saturated_events.values())
            ),
            n_pedestal=self.n_pedestal,
            n_survived_pedestals=self.n_survived_pedestals,
            n_triggered_tel1=n_triggered[0],
            n_triggered_tel2=n_triggered[1],
            swat_event_ids_used=getattr(source, "swat_event_ids_available", False),
        ))

        if self.wobble_in_output_name and wobble is not None and not pointing_manual:
            new_path = output_path.with_name(add_to_file_name(output_path.name, wobble))
            os.replace(output_path, new_path)
            self.log.info("Output file renamed to %s", new_path)


def add_to_file_name(file_name, text):
    """
    Add ``text`` to a file name, before its extensions,
    e.g. add_to_file_name("events.dl1.h5", "W1") == "events_W1.dl1.h5"
    """
    stem, dot, extensions = file_name.partition(".")
    return f"{stem}_{text}{dot}{extensions}"

def main():
    processor = ProcessorTool()
    processor.run()

if __name__ == '__main__':

    main()
