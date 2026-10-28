import os
from collections import Counter

import numpy as np
from tqdm import tqdm
from ctapipe.calib import CameraCalibrator
from ctapipe.containers import SchedulingBlockType
from ctapipe.core import Tool
from ctapipe.core.traits import Bool, ComponentName, List, UseEnum, classes_with_traits, flag
from ctapipe.image import ImageProcessor
from ctapipe.image.cleaning import NSBImageCleaner
from ctapipe.io import EventSource, DataWriter, MonitoringSource
from ctapipe.reco import ShowerProcessor

from sst1mpipe.calib import R0R1Calibrator, ImageSaturationCorrector
from sst1mpipe.io import attach_monitoring_file, write_dl1_info
from sst1mpipe.io.containers import PEDESTAL_EVENT_TYPES
from sst1mpipe.io.monitoring_source import SST1MMonitoringSource
from sst1mpipe.utils.cleaning import DBSCANImageCleaner, TimeDBSCANImageCleaner
from sst1mpipe.io.zmq_event_source import ZMQEventSource


class ProcessorTool(Tool):
    """
    Process data from lower-data levels up to DL1 including image
    extraction and optionally image parameterization.
    This implementation is based on the ctapipe.tool.ProcessorTool

    For the SST-1M raw data (R0, SST1MEventSource), the R0 -> R1 calibration is done
    by the R0R1Calibrator, with the statistics of the ADC samples of the pedestal events
    (voltage drop correction and dead pixels). The statistics of the calibrated images of
    the pedestal events are used by the NSBImageCleaner, which raises the picture threshold
    of the pixels with a high pedestal std.

    These statistics and the pointing of the telescope are given by the monitoring sources
    (monitoring_source_list), e.g. the monitoring file of the run written by
    sst1mpipe-monitoring (SST1MMonitoringSource), as in ctapipe-process: the monitoring
    containers of each event are filled before its calibration. The events without monitoring
    data (no chunk of pedestal events at their time) are skipped. The monitoring data are
    attached to the output file (attach_monitoring). Without monitoring source, the raw data are
    calibrated without these statistics, and the output file has no pointing.

    The charges and peak times of the saturated pixels are corrected by the ImageSaturationCorrector.
    How each telescope is calibrated is set in the configuration: for the simulated telescopes
    1 and 2, the R1 waveforms of the simulation are kept and corrected for the PDE drop
    (R0R1Calibrator.calibrate_r0 false, pde_drop_factor), without saturation correction
    (ImageSaturationCorrector.apply_correction false).
    Only the events of the runs with an allowed scheduling block type (allowed_sb_types)
    are processed, e.g. the observations but not the transitions between two wobbles.

    The pedestal events are used by the monitoring (sst1mpipe-monitoring) and are not written:
    they are calibrated and cleaned to count the fraction of them which survive the cleaning
    (n_survived_pedestals in /dl1/info). The average fraction of pixels whose picture threshold
    is raised by the NSBImageCleaner is also given in /dl1/info (raised_threshold_fraction).
    The shower geometry is reconstructed by the ShowerProcessor if the DL2 is written
    (DataWriter.write_dl2).
    """

    name = 'sst1mpipe-process'
    description = __doc__
    examples = """
    Simulations:

    > sst1mpipe-process -i mysim.simtel.gz -o events.dl1.h5 --config sst1mpipe_mc_config_low_nsb.json

    Raw data of a run with its monitoring file (sst1mpipe-monitoring):

    > sst1mpipe-process -i run.fits.fz -o run.dl1.h5 --config sst1mpipe_data_config.json
      --monitoring-source SST1MMonitoringSource --SST1MMonitoringSource.input_files=run.monitoring.h5

    ZMQ stream:

    > sst1mpipe-process -i tcp://localhost:24593 -o events.dl1.h5 --config sst1mpipe_rta_config.json
    """

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

    monitoring_source_list = List(
        ComponentName(MonitoringSource),
        default_value=[],
        help=(
            "Monitoring sources filling the monitoring containers of the events before their"
            " calibration, e.g. SST1MMonitoringSource (monitoring file of sst1mpipe-monitoring)."
            " Later sources overwrite earlier ones if their monitoring types overlap."
        ),
    ).tag(config=True)

    attach_monitoring = Bool(
        default_value=True,
        help=(
            "Copy the monitoring data of the monitoring files of the monitoring sources into the"
            " output file (pointing, statistics of the pedestal events)"
        ),
    ).tag(config=True)

    aliases = {
        ("i", "input"): "EventSource.input_url",
        ("o", "output"): "DataWriter.output_path",
        ("t", "allowed-tels"): "EventSource.allowed_tels",
        ("m", "max-events"): "EventSource.max_events",
        "monitoring-source": "ProcessorTool.monitoring_source_list",
    }

    flags = {

        **flag(
            "progress",
            "ProcessorTool.progress_bar",
            "show a progress bar during event processing",
            "don't show a progress bar during event processing",)
    }

    classes = (
        [
            CameraCalibrator, ImageProcessor, DataWriter, DBSCANImageCleaner, TimeDBSCANImageCleaner,
            R0R1Calibrator, ImageSaturationCorrector, ShowerProcessor,
            SST1MMonitoringSource,
        ]
        + classes_with_traits(EventSource)
        + classes_with_traits(MonitoringSource)
    )

    def setup(self):

        if ZMQEventSource.is_compatible(self.config.EventSource.input_url): # temporary fix since tcp:// url is not accepted by EventSource
            self.event_source = self.enter_context(ZMQEventSource(parent=self))
        else:
            self.event_source = self.enter_context(EventSource(parent=self))
        subarray = self.event_source.subarray
        # monitoring data of the events, e.g. the monitoring file of the run (sst1mpipe-monitoring)
        self.monitoring_sources = [
            self.enter_context(MonitoringSource.from_name(name, subarray=subarray, parent=self))
            for name in self.monitoring_source_list
        ]
        self._monitoring_pointing = any(getattr(source, "has_pointings", False) for source in self.monitoring_sources)
        # R0 -> R1 calibration of the raw data, PDE drop correction of the simulations: set per
        # telescope in the configuration (R0R1Calibrator.calibrate_r0)
        self.r0_r1_calibrator = R0R1Calibrator(parent=self, subarray=subarray)
        self.image_saturation_corrector = ImageSaturationCorrector(parent=self, subarray=subarray)
        self.camera_calibrator = CameraCalibrator(parent=self, subarray=self.event_source.subarray)
        self.image_processor = ImageProcessor(parent=self, subarray=self.event_source.subarray)
        # the writer is closed in finish(), to read back the output file. If the processing
        # fails before, it is closed when the tool exits.
        self.writer = DataWriter(event_source=self.event_source, parent=self)
        # reconstruction of the shower geometry, written in the DL2
        self.shower_processor = ShowerProcessor(parent=self, subarray=self.event_source.subarray)
        self._sb_types = {}
        self.n_skipped_events = 0
        self.n_events_without_monitoring = 0
        # processing summary, counted in the event loop
        self.n_triggered = Counter()
        self.n_pedestal = 0
        self.n_survived_pedestals = 0
        # sum and number of the fractions of pixels with a raised picture threshold
        self._raised_threshold = [0.0, 0]
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
            self.n_triggered.update(event.trigger.tels_with_trigger)
            is_pedestal = event.trigger.event_type in PEDESTAL_EVENT_TYPES
            self.n_pedestal += is_pedestal
            if len(self.monitoring_sources) > 0 and not self.fill_monitoring(event):
                self.n_events_without_monitoring += not is_pedestal
                continue

            # R0 -> R1 calibration with the statistics of the monitoring sources (raw data),
            # PDE drop correction of the R1 waveforms (simulations)
            self.r0_r1_calibrator(event)
            self.camera_calibrator(event)
            # the saturated pixels are corrected with the R0 waveforms
            self.image_saturation_corrector(event)
            self.image_processor(event)
            self.add_raised_threshold_fraction(event)
            # the pedestal events are not written: fraction of them which survive the cleaning
            if is_pedestal:
                self.n_survived_pedestals += any(
                    np.isfinite(dl1.parameters.hillas.intensity) for dl1 in event.dl1.tel.values()
                )
                continue
            if self.writer.write_dl2:
                self.shower_processor(event)
            self.writer(event)

    def fill_monitoring(self, event):
        """
        Fill the monitoring containers of the event with the monitoring sources.
        False if the event has no monitoring data at its time (statistics of the pedestal events,
        pointing if given by the monitoring sources): it cannot be calibrated.
        """
        for monitoring_source in self.monitoring_sources:
            monitoring_source.fill_monitoring_container(event)
        for tel_id in event.r0.tel:
            monitoring = event.monitoring.tel[tel_id]
            statistics = monitoring.camera.pixel_statistics
            no_pointing = self._monitoring_pointing and not np.isfinite(monitoring.pointing.altitude)
            if statistics.pedestal_waveform.std is None or statistics.pedestal_image.std is None or no_pointing:
                return False
        return True

    def add_raised_threshold_fraction(self, event):
        """
        Fraction of the pixels of the cleaned images whose picture threshold is raised by the
        NSBImageCleaner (pedestal_factor * std of the pedestal images above picture_threshold_pe)
        """
        cleaner = self.image_processor.clean
        if not isinstance(cleaner, NSBImageCleaner):
            return
        for tel_id in event.dl1.tel:
            std = event.monitoring.tel[tel_id].camera.pixel_statistics.pedestal_image.std
            if std is None:
                continue
            raised = cleaner.pedestal_factor.tel[tel_id] * np.asarray(std) > cleaner.picture_threshold_pe.tel[tel_id]
            self._raised_threshold[0] += np.mean(raised)
            self._raised_threshold[1] += 1

    @property
    def raised_threshold_fraction(self):
        """Average fraction of pixels with a raised picture threshold, None without statistics"""
        total, n = self._raised_threshold
        return total / n if n > 0 else None

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

    def finish(self):

        if self.event_source.is_simulation:
            self.writer.write_simulated_shower_distributions(self.event_source.simulated_shower_distributions)
        self._close_writer()
        output_path = self.writer.output_path
        if self.attach_monitoring:
            for monitoring_source in self.monitoring_sources:
                for monitoring_file in getattr(monitoring_source, "input_files", []):
                    attach_monitoring_file(monitoring_file, output_path)
                    self.log.info("Monitoring data of %s attached to the output file", monitoring_file)
        if self.n_skipped_events > 0:
            self.log.warning("%d events skipped (scheduling block type not allowed)", self.n_skipped_events)
        if self.n_events_without_monitoring > 0:
            self.log.warning(
                "%d events skipped (no monitoring data at their time)", self.n_events_without_monitoring,
            )

        n_triggered = {int(tel_id): self.n_triggered[tel_id] for tel_id in self.event_source.subarray.tel_ids}
        for tel_id, n in n_triggered.items():
            self.log.info("Number of triggered events of telescope %d: %d", tel_id, n)
        self.log.info("Number of pedestal events: %d", self.n_pedestal)
        if self.n_pedestal > 0:
            self.log.info(
                "Fraction of pedestal events that survived cleaning: %f", self.n_survived_pedestals / self.n_pedestal,
            )
        if self.raised_threshold_fraction is not None:
            self.log.info(
                "Average (per event) fraction of pixels with raised picture threshold: %f",
                self.raised_threshold_fraction,
            )

        # observation information, from the event source (SST1MEventSource)
        source = self.event_source
        target = getattr(source, "target", None)
        wobble = getattr(source, "wobble", None)
        pointing = getattr(source, "pointing", None)
        pointing_manual = getattr(source, "pointing_manual", False)
        self.log.info("Target: %s, wobble: %s, pointing: %s (manual: %s)", target, wobble, pointing, pointing_manual)

        # files of the telescopes calibrated from R0 (none for the simulations)
        calibration_files, window_files = None, None
        tel_ids = [tel_id for tel_id, n in n_triggered.items() if n > 0 and self.r0_r1_calibrator.calibrate_r0.tel[tel_id]]
        if len(tel_ids) > 0:
            calibration_files = ",".join(str(self.r0_r1_calibrator.calibration_file_path(t)) for t in tel_ids)
            window_files = ",".join(str(self.r0_r1_calibrator.window_transmittance_file_path(t)) for t in tel_ids)

        # telescopes 1 and 2: 21 and 22, or 1 and 2 in the simulations
        n_triggered_tel = {1: 0, 2: 0}
        for tel_id, n in n_triggered.items():
            n_triggered_tel[tel_id - 20 if tel_id > 20 else tel_id] = n
        write_dl1_info(output_path, dict(
            calib_file=calibration_files,
            window_file=window_files,
            target=target,
            wobble=wobble,
            ra=None if pointing is None else pointing.ra.deg,
            dec=None if pointing is None else pointing.dec.deg,
            manual_coords=pointing_manual,
            n_saturated=sum(self.image_saturation_corrector.n_saturated_events.values()),
            n_pedestal=self.n_pedestal,
            n_survived_pedestals=self.n_survived_pedestals,
            raised_threshold_fraction=self.raised_threshold_fraction,
            n_triggered_tel1=n_triggered_tel[1],
            n_triggered_tel2=n_triggered_tel[2],
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
