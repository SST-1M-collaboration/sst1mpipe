from tqdm import tqdm
from ctapipe.calib import CameraCalibrator
from ctapipe.core import Tool
from ctapipe.core.traits import Bool, flag
from ctapipe.image import ImageProcessor
from ctapipe.io import EventSource, DataWriter

from sst1mpipe.utils.cleaning import DBSCANImageCleaner, TimeDBSCANImageCleaner
from sst1mpipe.io.zmq_event_source import ZMQEventSource


class ProcessorTool(Tool):
    """
    Process data from lower-data levels up to DL1 including image
    extraction and optionally image parameterization.
    This implementation is based on the ctapipe.tool.ProcessorTool
    """

    name = 'sst1mpipe-process'
    description = __doc__
    examples = ("sst1mpipe-process -i mysim.simtel.gz -o events.dl1.h5",
                "sst1mpipe-process -i tcp://localhost:24593 -o events.dl1.h5 "
                "--config sst1mpipe/data/sst1mpipe_rta_config.json --log-level INFO")

    progress_bar = Bool(
        help="show progress bar during processing", default_value=False
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

    classes = [DBSCANImageCleaner, TimeDBSCANImageCleaner]

    def setup(self):

        if ZMQEventSource.is_compatible(self.config.EventSource.input_url):
            self.event_source = self.enter_context(ZMQEventSource(parent=self))
        else:
            self.event_source = self.enter_context(EventSource(parent=self))
        self.camera_calibrator = CameraCalibrator(parent=self, subarray=self.event_source.subarray)
        self.image_processor = ImageProcessor(parent=self, subarray=self.event_source.subarray)
        self.writer = self.enter_context(DataWriter(event_source=self.event_source, parent=self))

    def start(self):

        for event in tqdm(
            self.event_source,
            desc=self.event_source.__class__.__name__,
            total=self.event_source.max_events,
            disable=not self.progress_bar,
        ):
            self.camera_calibrator(event)
            self.image_processor(event)
            self.writer(event)

    def finish(self):

        pass

def main():
    processor = ProcessorTool()
    processor.run()

if __name__ == '__main__':

    main()
