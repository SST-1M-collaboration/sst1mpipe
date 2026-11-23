r"""
Datasets: the list of the runs of an analysis, in a text file (ECSV) easy to read in Python and in bash,
to reproduce an analysis and archive which data were used (with the configurations, e.g. in a zip).

One line per run and telescope: ``obs_id`` (see `sst1mpipe.io.sst1m_event_source.observation_ids`),
``tel_id`` and the raw file relative to the root of the data, given in the header with the computer
of the data (``host``, ``root``)::

    # %ECSV 1.0
    # ---
    # ...
    # meta: !!omap
    # - {name: crab_paper_2025}
    # - {host: calculus}
    # - {root: /net}
    # ...
    obs_id tel_id file
    2024011010012 21 cs1/data/raw/2024/01/10/SST1M1_20240110_0012.fits.fz

In bash::

    ROOT=$(sed -n 's/^# - {root: \(.*\)}$/\1/p' crab_paper_2025.ecsv)
    grep -v '^#' crab_paper_2025.ecsv | tail -n +2 | while read obs_id tel_id file; do echo "${ROOT}/${file}"; done
"""
import datetime
import logging
from pathlib import Path

import numpy as np
from astropy.table import Table

from sst1mpipe.io.sst1m_event_source import observation_ids, parse_file_name, read_target, tel_id_from_file_name

__all__ = ["dataset_paths", "find_runs", "make_dataset", "read_dataset", "write_dataset"]

logger = logging.getLogger(__name__)


def find_runs(root, night, telescopes=(1, 2)):
    """
    Raw files of a night (``YYYYMMDD``) of the telescopes: ``<root>/cs<n>/data/raw/YYYY/MM/DD/`` and
    its subdirectories, ``SST1M<n>_<night>_*.fits.fz``, sorted
    """
    night = str(night)
    files = []
    for telescope in telescopes:
        directory = Path(root) / f"cs{telescope}" / "data" / "raw" / night[:4] / night[4:6] / night[6:8]
        files += sorted(directory.rglob(f"SST1M{telescope}_*.fits.fz"))
    return files


def make_dataset(root, nights, target=None, telescopes=(1, 2), name="", description="", host="", log=None):
    """
    Dataset of the runs of the nights of ``root`` whose target (TARGET field of the raw file, without
    the wobble, compared in lower case) is ``target`` (all the runs if None)

    Returns
    -------
    astropy.table.Table: ``obs_id``, ``tel_id``, ``file`` (relative to ``root``), sorted by obs_id;
    ``meta``: name, description, host (computer of the data, e.g. calculus), root, target, nights,
    creation date and sst1mpipe version

    ``log``: logger of the files read (the logger of this module by default)
    """
    import sst1mpipe

    log = log or logger
    root = Path(root).absolute()
    rows = []
    for night in nights:
        paths = find_runs(root, night, telescopes)
        if not paths:
            night = str(night)
            log.warning("night %s: no raw files in %s", night, ", ".join(
                str(root / f"cs{t}" / "data" / "raw" / night[:4] / night[4:6] / night[6:8]) for t in telescopes))
        for path in paths:
            date_run, tel_id = parse_file_name(path), tel_id_from_file_name(path)
            if date_run is None or tel_id is None:
                log.warning("%s: not a raw file name SST1M<n>_<date>_<run>, skipped", path)
                continue
            try:
                run_target = read_target(path)
            except Exception as error:  # noqa: BLE001 (broken file)
                log.warning("%s: not readable (%s: %s), skipped", path, type(error).__name__, error)
                continue
            if target is not None and (run_target or "").lower() != target.lower():
                log.info("%s: target %s", path, run_target)
                continue
            _, obs_id = observation_ids(date_run[0], tel_id, date_run[1])
            log.info("%s: target %s, obs_id %d: in the dataset", path, run_target, obs_id)
            rows.append((obs_id, tel_id, str(path.relative_to(root))))
        log.info("night %s: %d runs in the dataset", night, sum(1 for row in rows if str(row[0]).startswith(str(night))))
    table = Table(rows=rows or None, names=["obs_id", "tel_id", "file"], dtype=[np.int64, np.int16, str])
    table.sort("obs_id")
    table.meta.update({
        "name": name, "description": description, "host": host, "root": str(root), "target": target or "",
        "nights": [int(night) for night in nights],
        "created": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
        "sst1mpipe_version": sst1mpipe.__version__,
    })
    return table


def write_dataset(table, path, overwrite=False):
    """Write a dataset (ECSV)"""
    table.write(path, format="ascii.ecsv", overwrite=overwrite)


def read_dataset(path):
    """
    Dataset of a file written by `write_dataset` (`make_dataset`), or of sst1mpipe by name
    (``sst1mpipe/resources/datasets/<name>.ecsv``, e.g. ``read_dataset("crab_paper_2025")``):
    ``meta["root"]`` is the root of the data
    """
    from sst1mpipe.resources import DATASETS_DIR

    path = Path(path)
    if not path.is_file() and (Path(DATASETS_DIR) / f"{path.name.removesuffix('.ecsv')}.ecsv").is_file():
        path = Path(DATASETS_DIR) / f"{path.name.removesuffix('.ecsv')}.ecsv"
    return Table.read(path, format="ascii.ecsv")


def dataset_paths(table, root=None):
    """Absolute paths of the files of a dataset, with its root or ``root`` (data moved)"""
    root = Path(root if root is not None else table.meta["root"])
    return [root / file for file in table["file"]]
