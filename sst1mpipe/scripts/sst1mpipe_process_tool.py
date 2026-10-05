import os

from tqdm import tqdm
from ctapipe.calib import CameraCalibrator
from ctapipe.containers import EventType
from ctapipe.core import Tool
from ctapipe.core.traits import Bool, flag
from ctapipe.image import ImageProcessor
from ctapipe.io import EventSource, DataWriter, SimTelEventSource

from sst1mpipe.calib import R0R1Calibrator
from sst1mpipe.io import compute_dl1_summary, get_used_qe_simtel, write_dl1_info
from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.utils.monitoring_pedestals import R0PedestalMonitor
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
    For the simulations (SimTelEventSource), the R1 waveforms are corrected for the
    PDE drop by the R0R1Calibrator (mc_pde_correction).
    """

    name = 'sst1mpipe-process'
    description = __doc__
    examples = ("sst1mpipe-process -i mysim.simtel.gz -o events.dl1.h5",
                "sst1mpipe-process -i tcp://localhost:24593 -o events.dl1.h5 "
                "--config sst1mpipe/data/sst1mpipe_rta_config.json --log-level INFO")

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

    classes = [DBSCANImageCleaner, TimeDBSCANImageCleaner, R0R1Calibrator, R0PedestalMonitor]

    def setup(self):

        if ZMQEventSource.is_compatible(self.config.EventSource.input_url): # temporary fix since tcp:// url is not accepted by EventSource
            self.event_source = self.enter_context(ZMQEventSource(parent=self))
        else:
            self.event_source = self.enter_context(EventSource(parent=self))
        # R0 -> R1 calibration of the SST-1M raw data, PDE drop correction of the simulations.
        # The other sources (ZMQ) provide calibrated R1 data
        self.r0_pedestal_monitor = None
        self.r0_r1_calibrator = None
        subarray = self.event_source.subarray
        if isinstance(self.event_source, SST1MEventSource):
            self.r0_pedestal_monitor = R0PedestalMonitor(parent=self, subarray=subarray)
            self.r0_r1_calibrator = R0R1Calibrator(parent=self, subarray=subarray)
        elif isinstance(self.event_source, SimTelEventSource):
            self.r0_r1_calibrator = R0R1Calibrator(
                parent=self, subarray=subarray, simulated_pde_files=get_used_qe_simtel(self.event_source),
            )
        self.camera_calibrator = CameraCalibrator(parent=self, subarray=self.event_source.subarray)
        self.image_processor = ImageProcessor(parent=self, subarray=self.event_source.subarray)
        # the writer is closed in finish(), to read back the output file. If the processing
        # fails before, it is closed when the tool exits.
        self.writer = DataWriter(event_source=self.event_source, parent=self)
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
            if self.r0_r1_calibrator is not None:
                self.calibrate_r0_r1(event)
            self.camera_calibrator(event)
            self.image_processor(event)
            self.writer(event)

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

    def finish(self):

        self._close_writer()
        output_path = self.writer.output_path

        # processing summary, computed from the content of the output file
        summary = compute_dl1_summary(output_path)
        self.log.info("Number of events: %d", summary["n_events"])
        for tel_id, n in summary["n_triggered"].items():
            self.log.info("Number of triggered events of telescope %d: %d", tel_id, n)
        self.log.info("Number of pedestal events: %d", summary["n_pedestal"])
        if summary["n_pedestal"] > 0:
            self.log.info(
                "Fraction of pedestal events that survived cleaning: %f",
                summary["n_survived_pedestals"] / summary["n_pedestal"],
            )

        # observation information, from the event source (SST1MEventSource)
        source = self.event_source
        target = getattr(source, "target", None)
        wobble = getattr(source, "wobble", None)
        pointing = getattr(source, "pointing", None)
        pointing_manual = getattr(source, "pointing_manual", False)
        self.log.info("Target: %s, wobble: %s, pointing: %s (manual: %s)", target, wobble, pointing, pointing_manual)

        calibration_files = None
        if self.r0_pedestal_monitor is not None:  # observed data
            calibration_files = ",".join(
                str(self.r0_r1_calibrator.calibration_file_path(tel_id))
                for tel_id, n in summary["n_triggered"].items() if n > 0
            )

        n_triggered = list(summary["n_triggered"].values()) + [0, 0]
        write_dl1_info(output_path, dict(
            calib_file=calibration_files,
            target=target,
            wobble=wobble,
            ra=None if pointing is None else pointing.ra.deg,
            dec=None if pointing is None else pointing.dec.deg,
            manual_coords=pointing_manual,
            n_saturated=None,  # the saturation correction is not applied yet
            n_pedestal=summary["n_pedestal"],
            n_survived_pedestals=summary["n_survived_pedestals"],
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
