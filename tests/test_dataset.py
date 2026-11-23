"""Datasets: list of the runs of an analysis, read in Python and in bash"""
import logging
import subprocess

from astropy.io import fits

from conftest import DARK_FILE, OBSERVATION_RUNS, TRANSITION_FILES
from sst1mpipe.io.dataset import dataset_paths, find_runs, read_dataset
from sst1mpipe.scripts.sst1mpipe_make_dataset import main as make_dataset


def raw_data_root(tmp_path):
    """Root of raw data with the test files: <root>/cs<n>/data/raw/YYYY/MM/DD/(SST1M<n>/)<file>"""
    root = tmp_path / "raw_data"
    files = [OBSERVATION_RUNS[21]["path"], OBSERVATION_RUNS[22]["path"], TRANSITION_FILES[21], DARK_FILE]
    for path in files:
        telescope, night = path.name[5], path.name.split("_")[1]
        directory = root / f"cs{telescope}" / "data" / "raw" / night[:4] / night[4:6] / night[6:8]
        # the files of tel 2 in a subdirectory SST1M2, as on some disks
        directory = directory / "SST1M2" if telescope == "2" else directory
        directory.mkdir(parents=True, exist_ok=True)
        (directory / path.name).symlink_to(path)
    return root


def test_find_runs(tmp_path):
    root = raw_data_root(tmp_path)
    assert [path.name for path in find_runs(root, 20260120)] == [
        "SST1M1_20260120_1179.fits.fz", "SST1M2_20260120_1102.fits.fz"]
    assert [path.name for path in find_runs(root, "20260120", telescopes=[2])] == ["SST1M2_20260120_1102.fits.fz"]
    assert find_runs(root, 20250101) == []


def test_make_dataset(tmp_path, capsys):
    """the runs of a target (not the transitions and dark runs), readable in Python and bash"""
    root = raw_data_root(tmp_path)
    output = tmp_path / "mrk421.ecsv"
    make_dataset(["--host", "calculus", "--root", str(root), "--target", "mrk421", "-o", str(output),
                  "--description", "test",
                  "20260119", "20260120", "20260121"])
    assert "2 runs" in capsys.readouterr().out

    dataset = read_dataset(output)
    assert dataset["obs_id"].tolist() == [2026012011179, 2026012021102]
    assert dataset["tel_id"].tolist() == [21, 22]
    assert dataset.meta["root"] == str(root) and dataset.meta["target"] == "mrk421"
    assert dataset.meta["host"] == "calculus"
    assert dataset.meta["nights"] == [20260119, 20260120, 20260121]
    assert [path.resolve() for path in dataset_paths(dataset)] == [
        OBSERVATION_RUNS[21]["path"].resolve(), OBSERVATION_RUNS[22]["path"].resolve()]
    assert dataset_paths(dataset, root="/moved")[0].parts[:2] == ("/", "moved")

    # bash: root from the header, one line per run
    script = (
        f"ROOT=$(sed -n 's/^# - {{root: \\(.*\\)}}$/\\1/p' {output}); "
        f"grep -v '^#' {output} | tail -n +2 | while read obs_id tel_id file; do echo \"$obs_id ${{ROOT}}/$file\"; done"
    )
    lines = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True).stdout.split("\n")
    assert lines[:2] == [f"{obs_id} {path}" for obs_id, path in zip(dataset["obs_id"], dataset_paths(dataset), strict=True)]

    # all the runs without target
    make_dataset(["--root", str(root), "-o", str(output), "--overwrite", "20260119", "20260120", "20260121"])
    assert len(read_dataset(output)) == 4


def test_old_and_broken_raw_files(tmp_path, caplog):
    """TARGET in a table of another name (old raw files), unreadable file skipped, verbose log"""
    from sst1mpipe.io.dataset import make_dataset as make

    root = raw_data_root(tmp_path)
    night = root / "cs1" / "data" / "raw" / "2026" / "01" / "20"
    old = fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(name="DataStream"),
                        fits.BinTableHDU(name="OldEvents", header=fits.Header({"TARGET": "Crab_W2_83.63_22.01"}))])
    old.writeto(night / "SST1M1_20260120_0001.fits.fz")
    (night / "SST1M1_20260120_0002.fits.fz").write_bytes(b"not a fits file")

    with caplog.at_level(logging.INFO, logger="sst1mpipe.io.dataset"):
        dataset = make(root, ["20260120"], target="crab")
    assert dataset["obs_id"].tolist() == [2026012010001]
    messages = caplog.text
    assert "SST1M1_20260120_0001.fits.fz: target Crab, obs_id 2026012010001: in the dataset" in messages
    assert "SST1M1_20260120_1179.fits.fz: target MRK421" in messages
    assert "SST1M1_20260120_0002.fits.fz: not readable" in messages
    assert "night 20260120: 1 runs in the dataset" in messages
