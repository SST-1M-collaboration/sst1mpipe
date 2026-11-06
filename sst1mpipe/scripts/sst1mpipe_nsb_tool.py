"""
NSB level of SST-1M runs and voltage drop of the SiPMs, from the statistics of the ADC samples of
the pedestal events of monitoring files (sst1mpipe-monitoring) or DL1 files (sst1mpipe-process).

For each chunk of pedestal events and telescope, from the variance of the ADC samples (mean over
the pixels of their std, squared), with the SiPMVoltageDropModel: baseline shift, gain drop,
charge drop (voltage drop correction of the R0R1Calibrator), NSB rate and NSB photon rate per
pixel. The table is written in an ECSV, FITS or HDF5 file (by the extension of the output file),
and optionally plotted vs. time.
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import QTable, vstack
from ctapipe.core import Provenance, Tool, ToolConfigurationError
from ctapipe.core.traits import Bool, List, Path, flag
from ctapipe.io import write_table

from sst1mpipe.utils.nsb_tools import (
    SiPMVoltageDropModel,
    mean_pedestal_std,
    pedestal_time,
    plot_nsb_rate,
    read_pedestal_tables,
)
from sst1mpipe.utils.utils import get_subarray

HDF5_TABLE = "/nsb"


class NSBTool(Tool):
    """NSB level and voltage drop of the SiPMs of SST-1M runs, from their pedestal events"""

    name = "sst1mpipe-nsb"
    description = __doc__
    examples = (
        "sst1mpipe-nsb -i SST1M1_20260120_1179.monitoring.h5 -o nsb.ecsv --plot nsb.png\n"
        "sst1mpipe-nsb -i run_1.dl1.h5 -i run_2.dl1.h5 -o nsb.h5 --SiPMVoltageDropModel.photon_detection_efficiency=0.3"
    )

    input_files = List(
        Path(exists=True, directory_ok=False),
        default_value=[],
        help="Monitoring files (sst1mpipe-monitoring) or DL1 files (sst1mpipe-process), one -i per file",
    ).tag(config=True)

    output_path = Path(
        default_value=None, allow_none=True, directory_ok=False,
        help=f"Output table: ECSV (.ecsv), FITS (.fits) or HDF5 (.h5, table {HDF5_TABLE})",
    ).tag(config=True)

    plot_path = Path(
        default_value=None, allow_none=True, directory_ok=False,
        help="Plot of the NSB rate vs. time of the telescopes (e.g. .png, .pdf), not plotted if not given",
    ).tag(config=True)

    overwrite = Bool(default_value=False, help="Overwrite the output files if they exist").tag(config=True)

    aliases = {
        ("i", "input"): "NSBTool.input_files",
        ("o", "output"): "NSBTool.output_path",
        ("p", "plot"): "NSBTool.plot_path",
    }

    flags = {
        **flag("overwrite", "NSBTool.overwrite", "overwrite the output files", "do not overwrite the output files"),
    }

    classes = [SiPMVoltageDropModel]

    def setup(self):
        if len(self.input_files) == 0:
            raise ToolConfigurationError("No input file (-i/--input or NSBTool.input_files)")
        if self.output_path is None:
            raise ToolConfigurationError("The output file is required (-o/--output or NSBTool.output_path)")
        if self.output_path.suffix not in (".ecsv", ".fits", ".h5", ".hdf5"):
            raise ToolConfigurationError(f"Unknown format of the output file {self.output_path}: .ecsv, .fits or .h5 expected")
        for path in (self.output_path, self.plot_path):
            if path is not None and path.exists() and not self.overwrite:
                raise ToolConfigurationError(f"Output file {path} exists, use --overwrite to overwrite it")
        self.model = SiPMVoltageDropModel(subarray=get_subarray(), parent=self)
        self.table = None

    def start(self):
        pedestals = read_pedestal_tables(self.input_files)
        if pedestals is None:
            raise ToolConfigurationError("No statistics of the pedestal events in the input files")
        for path in self.input_files:
            Provenance().add_input_file(str(path), role="Monitoring")

        rows = []
        for group in pedestals.group_by("tel_id").groups:
            tel_id = int(group["tel_id"][0])
            variance = mean_pedestal_std(group) ** 2
            n_events = np.asarray(group["n_events"])
            rows.append(QTable(dict(
                tel_id=group["tel_id"],
                time=pedestal_time(group),
                time_start=group["time_start"],
                time_end=group["time_end"],
                # pedestal events of the chunk (maximum over the pixels)
                n_events=n_events.reshape(len(n_events), -1).max(axis=1),
                pedestal_variance=variance,
                baseline_shift=self.model.baseline_shift(tel_id, variance),
                gain_drop=self.model.gain_drop(tel_id, variance),
                charge_drop=self.model.charge_drop(tel_id, variance),
                nsb_rate=self.model.nsb_rate(tel_id, variance),
                nsb_photon_rate=self.model.nsb_photon_rate(tel_id, variance),
            )))
            self.log.info("Telescope %d: %d chunks of pedestal events", tel_id, len(group))
        self.table = vstack(rows)
        self.pedestals = pedestals

    def finish(self):
        if self.output_path.suffix in (".h5", ".hdf5"):
            write_table(self.table, self.output_path, HDF5_TABLE, overwrite=self.overwrite)
        else:
            self.table.write(self.output_path, overwrite=self.overwrite)
        Provenance().add_output_file(str(self.output_path), role="NSB")
        self.log.info("NSB table written in %s", self.output_path)

        if self.plot_path is not None:
            fig, ax = plt.subplots(figsize=(10, 5))
            for group in self.pedestals.group_by("tel_id").groups:
                tel_id = int(group["tel_id"][0])
                plot_nsb_rate(group, tel_id, model=self.model, ax=ax, label=f"tel {tel_id}")
            ax.legend()
            fig.savefig(self.plot_path)
            plt.close(fig)
            self.log.info("NSB plot written in %s", self.plot_path)


def main():
    tool = NSBTool()
    tool.run()


if __name__ == "__main__":
    main()
