r"""
``sst1mpipe-make-dataset``: dataset of an analysis (`sst1mpipe.io.dataset`), the runs of some nights
of one target, in a text file (ECSV) with the computer and the root of the data:

    sst1mpipe-make-dataset --host calculus --root /net --target Crab --name crab_paper_2025 \
        -o crab_paper_2025.ecsv 20231017 20231108 20231218

Read in bash:

    ROOT=$(sed -n 's/^# - {root: \(.*\)}$/\1/p' crab_paper_2025.ecsv)
    grep -v '^#' crab_paper_2025.ecsv | tail -n +2 | while read obs_id tel_id file; do echo "${ROOT}/${file}"; done
"""
import argparse
from pathlib import Path

from sst1mpipe.io.dataset import make_dataset, write_dataset


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("nights", nargs="+", help="nights YYYYMMDD (directories of the raw files)")
    parser.add_argument("--root", required=True, type=Path,
                        help="root of the raw data: <root>/cs<n>/data/raw/YYYY/MM/DD/")
    parser.add_argument("--host", default="", help="computer of the data, e.g. calculus (written in the dataset)")
    parser.add_argument("--target", help="target of the runs (TARGET field of the raw files, without the wobble)")
    parser.add_argument("--telescopes", type=int, nargs="+", default=[1, 2], help="telescopes (1, 2)")
    parser.add_argument("--name", default="", help="name of the dataset")
    parser.add_argument("--description", default="", help="description of the dataset")
    parser.add_argument("-o", "--output", required=True, type=Path, help="dataset file (.ecsv)")
    parser.add_argument("--overwrite", action="store_true")
    options = parser.parse_args(args)

    name = options.name or options.output.stem
    table = make_dataset(options.root, options.nights, target=options.target, telescopes=options.telescopes,
                         name=name, description=options.description, host=options.host)
    write_dataset(table, options.output, overwrite=options.overwrite)
    for tel_id in sorted(set(table["tel_id"])):
        print(f"tel {tel_id}: {sum(table['tel_id'] == tel_id)} runs")
    print(f"{len(table)} runs in {options.output}")


if __name__ == "__main__":
    main()
