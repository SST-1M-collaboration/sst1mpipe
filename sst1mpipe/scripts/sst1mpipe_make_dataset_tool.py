r"""
``sst1mpipe-make-dataset``: dataset of an analysis (`sst1mpipe.io.dataset`), the runs of some nights
of one target, in a text file (ECSV) with the computer and the root of the data:

    sst1mpipe-make-dataset --host calculus --root /net --target Crab --name crab_paper_2025 \
        -o crab_paper_2025.ecsv 20231017 20231108 20231218

The nights are positional arguments (or ``MakeDataset.nights`` in a configuration file).

Read in bash:

    ROOT=$(sed -n 's/^# - {root: \(.*\)}$/\1/p' crab_paper_2025.ecsv)
    grep -v '^#' crab_paper_2025.ecsv | tail -n +2 | while read obs_id tel_id file; do echo "${ROOT}/${file}"; done
"""
import logging

from ctapipe.core import Tool, ToolConfigurationError
from ctapipe.core.traits import Bool, Int, List, Path, Unicode, flag
from traitlets import default

from sst1mpipe.io.dataset import make_dataset, write_dataset

__all__ = ["MakeDataset"]


class MakeDataset(Tool):
    """Dataset of the runs of some nights of one target (one line per run and telescope)"""

    name = "sst1mpipe-make-dataset"
    description = __doc__
    examples = (
        "sst1mpipe-make-dataset --host calculus --root /net --target Crab -o crab.ecsv 20240110 20240111"
    )

    root = Path(
        exists=True, file_ok=False, allow_none=True, default_value=None,
        help="Root of the raw data: <root>/cs<n>/data/raw/YYYY/MM/DD/",
    ).tag(config=True)
    nights = List(
        Unicode(), default_value=[],
        help="Nights YYYYMMDD (directories of the raw files), also given as positional arguments",
    ).tag(config=True)
    target = Unicode(
        default_value=None, allow_none=True,
        help="Target of the runs (TARGET field of the raw files, without the wobble, any case); all the runs if None",
    ).tag(config=True)
    telescopes = List(Int(), default_value=[1, 2], help="Telescopes (1, 2)").tag(config=True)
    dataset_name = Unicode("", help="Name of the dataset (default: name of the output file)").tag(config=True)
    dataset_description = Unicode("", help="Description of the dataset").tag(config=True)
    host = Unicode("", help="Computer of the data, e.g. calculus (written in the dataset)").tag(config=True)
    output_path = Path(
        directory_ok=False, allow_none=True, default_value=None, help="Dataset file (.ecsv)",
    ).tag(config=True)
    overwrite = Bool(False, help="Overwrite the output file").tag(config=True)

    aliases = {
        "root": "MakeDataset.root",
        "target": "MakeDataset.target",
        "telescopes": "MakeDataset.telescopes",
        "name": "MakeDataset.dataset_name",
        "description": "MakeDataset.dataset_description",
        "host": "MakeDataset.host",
        ("o", "output"): "MakeDataset.output_path",
    }
    flags = {
        **flag("overwrite", "MakeDataset.overwrite", "Overwrite the output file", "Do not overwrite the output file"),
    }

    @default("log_level")
    def _default_log_level(self):
        """verbose: the files read and the runs selected (--log-level WARNING: the warnings only)"""
        return logging.INFO

    def setup(self):
        # nights as positional arguments, as the input files of ctapipe-merge
        self.nights = [*self.nights, *self.extra_args]
        if self.root is None or self.output_path is None or not self.nights:
            raise ToolConfigurationError("--root, --output and nights (positional arguments) are required")
        invalid = [night for night in self.nights if not (len(night) == 8 and night.isdigit())]
        if invalid:
            raise ToolConfigurationError(f"Nights must be YYYYMMDD: {invalid}")
        if self.output_path.exists() and not self.overwrite:
            raise ToolConfigurationError(f"{self.output_path} exists: --overwrite to replace it")

    def start(self):
        self.dataset = make_dataset(
            self.root, self.nights, target=self.target, telescopes=self.telescopes,
            name=self.dataset_name or self.output_path.stem, description=self.dataset_description, host=self.host,
            log=self.log,
        )

    def finish(self):
        write_dataset(self.dataset, self.output_path, overwrite=self.overwrite)
        for tel_id in sorted(set(self.dataset["tel_id"])):
            self.log.info("tel %d: %d runs", tel_id, sum(self.dataset["tel_id"] == tel_id))
        self.log.info("%d runs in %s", len(self.dataset), self.output_path)


def main():
    MakeDataset().run()
