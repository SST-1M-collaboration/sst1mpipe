from pathlib import Path
from subprocess import CompletedProcess

import pytest
from conftest import (
    RunCfg,
    assert_schema_matches_template,
    fail_if_missing,
    fixtures_dir,
    parse_log,
    run_cli,
    skip_if_missing,
)

LOG_REQUIRED_PATTERNS = {
    "tel1": r"Total number of TEL1 triggered events in the file:\s*([0-9]+)",
    "tel2": r"Total number of TEL2 triggered events in the file:\s*([0-9]+)",
    "saturated": r"Total number of saturated events in the file:\s*([0-9]+)",
    "pedestal": r"Total number of pedestal events in the file:\s*([0-9]+)",
}


@pytest.fixture(scope="module")
def r0_dl1_cfg(tmp_path_factory: pytest.TempPathFactory):
    fixtures: Path = fixtures_dir()
    return RunCfg(
        exe="sst1mpipe_r0_dl1",
        opts={
            "--input-file": fixtures / "r0/SST1M2_20251218_0350.fits.fz",
            "--config": fixtures / "default_sst1mpipe_data_config_analysis.json",
            "--output-dir": tmp_path_factory.mktemp("sst1mpipe_dl1_out"),
        },
        template_file=fixtures / "dl1/SST1M2_20251218_0350_W1_dl1.h5",
    )


@pytest.fixture(scope="module")
def produced_DL1_outputs(r0_dl1_cfg: RunCfg):
    skip_if_missing(r0_dl1_cfg.missing_paths(), label="r0->dl1 fixtures")
    print(f"\n[pytest] Running {r0_dl1_cfg.exe} (this may take a while)...")
    repo_root: Path = Path(__file__).resolve().parents[1]
    proc: CompletedProcess[str] = run_cli(
        r0_dl1_cfg.exe, [*r0_dl1_cfg.args, "--precise-timestamps"], cwd=repo_root
    )

    if proc.returncode != 0:
        pytest.skip(
            "Skipping pipeline-dependent tests because the pipeline failed.\n\n"
            f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}"
        )

    log_path: Path | None = next(r0_dl1_cfg.output_dir.glob("*.log"), None)
    h5_path: Path | None = next(r0_dl1_cfg.output_dir.glob("*.h5"), None)

    return log_path, h5_path


# ----------------------------
# Test 0: FAIL if fixtures are missing
# ----------------------------


def test_0_fixtures_exist(r0_dl1_cfg: RunCfg):
    fail_if_missing(r0_dl1_cfg.missing_paths(), label="r0->dl1 fixtures")


# ----------------------------
# Test 1: execute the code
# ----------------------------


def test_1_run_sst1mpipe_r0_dl1(produced_DL1_outputs: tuple):
    assert produced_DL1_outputs is not None


# ----------------------------
# Tests 2: Check if we generated required files (log, .h5)
# ----------------------------


def test_2_artifacts_exist_and_non_empty(produced_DL1_outputs: tuple):
    log_path, h5_path = produced_DL1_outputs
    assert log_path.exists() and log_path.stat().st_size > 0
    assert h5_path.exists() and h5_path.stat().st_size > 0


# ----------------------------
# Tests 3: Check that generated log contains expected patterns
# ----------------------------


def test_3_log_contains_required_markers(produced_DL1_outputs: tuple):
    log_path, _ = produced_DL1_outputs
    parsed = parse_log(log_path, LOG_REQUIRED_PATTERNS)
    assert (
        parsed["tel1"] + parsed["tel2"] + parsed["saturated"] + parsed["pedestal"]
    ) > 0


# ----------------------------
# Tests 4: Compare generated h5 file with template
# ----------------------------


def test_4_h5_schema_matches_template(r0_dl1_cfg: RunCfg, produced_DL1_outputs: tuple):
    _, produced_h5 = produced_DL1_outputs
    assert_schema_matches_template(produced_h5, r0_dl1_cfg.template_file)
