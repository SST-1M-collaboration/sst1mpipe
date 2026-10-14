from collections import Counter

import numpy as np
from ctapipe.core import Container, Field, Tool
from ctapipe.core.traits import Bool, Int, Path, flag
from ctapipe.io import EventSource, HDF5TableWriter
from tqdm import tqdm

from sst1mpipe.io.sst1m_event_source import SST1MEventSource
from sst1mpipe.trigger import TriggerPatch7Emulator, TriggerTDSCANEmulator
from sst1mpipe.utils import add_event_id


class TriggerEmulationContainer(Container):
    """Trigger emulation results of one telescope event."""

    obs_id = Field(-1, "Observation id")
    event_id = Field(-1, "Event id")
    tel_id = Field(-1, "Telescope id")
    event_type = Field(-1, "Camera event type (real data, 8 = pedestal), -1 for simulations")
    trigger_traces = Field(None, "Triplet traces, the FADC output (432, n_samples), patch_sw_id order")
    patch7_max = Field(-1, "Largest 7-patch cluster sum")
    patch7 = Field(False, "patch7 fired")
    tdscan_max = Field(np.nan, "Largest TDSCAN score")
    tdscan = Field(False, "TDSCAN fired")
    patch7_output = Field(None, "patch7 output per patch and sample (432, n_samples)")
    tdscan_output = Field(None, "TDSCAN output per patch and sample (432, n_samples)")
    camera_patch7 = Field(False, "patch7 fired in the camera (recorded output)")
    camera_patch7_mismatches = Field(-1, "Patches x samples where the emulated and the recorded patch7 outputs differ")


class TriggerTool(Tool):
    """
    Emulate the SST-1M camera trigger (patch7 and TDSCAN) on R0 data, real
    or simulated, and write the trigger traces and the trigger decisions.

    Real data can also be compared with the patch7 output recorded by the
    camera (--compare-camera).
    """

    name = "sst1mpipe-trigger"
    description = __doc__
    examples = (
        "sst1mpipe-trigger -i gamma.simtel.gz -o gamma_trigger.h5 --config sst1mpipe/data/sst1mpipe_trigger_emulator_mc.json\n"
        "sst1mpipe-trigger -i SST1M1_20260426_0250.fits.fz -o trigger.h5 "
        "--config sst1mpipe/data/sst1mpipe_trigger_emulator_data.json --compare-camera"
    )

    input_url = Path(exists=True, directory_ok=False, help="Input file: simtel or SST-1M .fits.fz").tag(config=True)
    output_path = Path(directory_ok=False, help="Output HDF5 file").tag(config=True)
    max_events = Int(None, allow_none=True, help="Maximum number of events to read").tag(config=True)
    overwrite = Bool(False, help="Overwrite the output file if it exists").tag(config=True)
    progress_bar = Bool(False, help="Show a progress bar during processing").tag(config=True)
    save_outputs = Bool(
        False, help="Also write the trigger outputs per patch and sample (432 x n_samples booleans per trigger)"
    ).tag(config=True)
    compare_camera = Bool(
        False, help="Real data: compare the emulated patch7 output with the one recorded by the camera"
    ).tag(config=True)

    aliases = {
        ("i", "input"): "TriggerTool.input_url",
        ("o", "output"): "TriggerTool.output_path",
        ("m", "max-events"): "TriggerTool.max_events",
    }

    flags = {
        **flag("progress", "TriggerTool.progress_bar",
               "show a progress bar during event processing", "don't show a progress bar"),
        **flag("overwrite", "TriggerTool.overwrite",
               "overwrite the output file", "don't overwrite the output file"),
        **flag("save-outputs", "TriggerTool.save_outputs",
               "write the trigger outputs per patch and sample", "write only the trigger decisions"),
        **flag("compare-camera", "TriggerTool.compare_camera",
               "compare with the patch7 output recorded by the camera", "don't compare with the camera"),
    }

    classes = [TriggerPatch7Emulator, TriggerTDSCANEmulator]

    def setup(self):
        if self.output_path.exists() and not self.overwrite:
            raise ValueError(f"{self.output_path} exists, use --overwrite to replace it")

        self.is_real_data = str(self.input_url).endswith(".fits.fz")
        if self.compare_camera and not self.is_real_data:
            raise ValueError("--compare-camera needs real data (.fits.fz): simulations have no recorded trigger output")
        if self.is_real_data:
            self.event_source = self.enter_context(
                SST1MEventSource([str(self.input_url)], max_events=self.max_events, parent=self))
        else:
            self.event_source = self.enter_context(
                EventSource(self.input_url, max_events=self.max_events, parent=self))

        subarray = self.event_source.subarray
        self.patch7 = TriggerPatch7Emulator(subarray=subarray, parent=self)
        self.tdscan = TriggerTDSCANEmulator(subarray=subarray, parent=self)

        self.writer = self.enter_context(HDF5TableWriter(self.output_path, group_name="trigger_emulation", mode="w"))
        for tel_id in subarray.tel_ids:
            table = f"tel_{tel_id:03d}"
            if not self.save_outputs:
                self.writer.exclude(table, "patch7_output")
                self.writer.exclude(table, "tdscan_output")
            if not self.compare_camera:
                self.writer.exclude(table, "camera_patch7")
                self.writer.exclude(table, "camera_patch7_mismatches")
        self.counts = Counter()

    def start(self):
        for event_number, event in enumerate(tqdm(
            self.event_source,
            desc=self.event_source.__class__.__name__,
            total=self.event_source.max_events,
            disable=not self.progress_bar,
        )):
            if self.is_real_data:
                # obs_id from the file name (SST1M1_<date>_<run>), as sst1mpipe_r0_dl1 does
                run_name = self.input_url.name.removesuffix(".fits.fz")
                event = add_event_id(event, filename=run_name, event_number=event_number)
            for tel_id in event.r0.tel.keys():
                self.writer.write(f"tel_{tel_id:03d}", [self.emulate(event, tel_id)])

    def emulate(self, event, tel_id):
        traces = self.patch7.trigger_traces(event, tel_id)
        patch7_output, patch7_scores = self.patch7.trigger_output(traces, tel_id)
        tdscan_output, tdscan_scores = self.tdscan.trigger_output(traces, tel_id)

        result = TriggerEmulationContainer(
            obs_id=event.index.obs_id,
            event_id=event.index.event_id,
            tel_id=tel_id,
            event_type=int(event.sst1m.r0.tel[tel_id].camera_event_type.value) if self.is_real_data else -1,
            trigger_traces=traces.astype(np.uint8),
            patch7_max=int(patch7_scores.max()),
            patch7=bool(patch7_output.any()),
            tdscan_max=float(tdscan_scores.max()),
            tdscan=bool(tdscan_output.any()),
            patch7_output=patch7_output,
            tdscan_output=tdscan_output,
        )
        self.counts[tel_id, "events"] += 1
        self.counts[tel_id, "patch7"] += result.patch7
        self.counts[tel_id, "tdscan"] += result.tdscan

        if self.compare_camera:
            camera_output = self.patch7.hardware_output(event, tel_id)
            result.camera_patch7 = bool(camera_output.any())
            result.camera_patch7_mismatches = int(np.sum(camera_output != patch7_output))
            self.counts[tel_id, "camera_patch7"] += result.camera_patch7
            self.counts[tel_id, "events_differing"] += result.camera_patch7_mismatches > 0
        return result

    def finish(self):
        for tel_id in sorted({tel for tel, _ in self.counts}):
            line = (f"tel {tel_id}: {self.counts[tel_id, 'events']} events, patch7 fired {self.counts[tel_id, 'patch7']}, "
                    f"TDSCAN fired {self.counts[tel_id, 'tdscan']}")
            if self.compare_camera:
                line += (f", camera patch7 fired {self.counts[tel_id, 'camera_patch7']}, "
                         f"events where emulated and camera patch7 outputs differ: {self.counts[tel_id, 'events_differing']}")
            print(line)


def main():
    tool = TriggerTool()
    tool.run()


if __name__ == "__main__":
    main()
