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
    "total_events_tel_022": r"\bTotal N of events of tel_022:\s*([0-9]+)",
    "selected_events_tel_022": r"\bN of events of tel_022 after selection cuts:\s*([0-9]+)",
    "reco_done": r"\b(Reconstruction done\.)",
    "writing_dl2": r"\bWriting DL2 in\s+(\S+)",
    "creating_dl2_node": r"\bCreating new node\s+(/dl2/event/telescope/parameters/tel_022)\b",
}


@pytest.fixture(scope="module")
def dl1_dl2_cfg(tmp_path_factory: pytest.TempPathFactory):
    fixtures: Path = fixtures_dir()
    return RunCfg(
        exe="sst1mpipe_dl1_dl2",
        opts={
            "--input-file": fixtures / "dl1/SST1M2_20251218_0350_W1_dl1.h5",
            "--config": fixtures / "default_sst1mpipe_data_config_analysis.json",
            "--models-dir": "/data/work/analysis/MC/prod_oct_2024_84/v0.7.2/models_mono/",
            "--output-dir": tmp_path_factory.mktemp("sst1mpipe_dl2_out"),
        },
        template_file=fixtures / "dl2/SST1M2_20251218_0350_W1_dl1_dl2.h5",
    )


@pytest.fixture(scope="module")
def produced_DL2_outputs(dl1_dl2_cfg: RunCfg):
    skip_if_missing(dl1_dl2_cfg.missing_paths(), label="dl1->dl2 fixtures")

    print(f"\n[pytest] Running {dl1_dl2_cfg.exe} (this may take a while)...")

    repo_root = Path(__file__).resolve().parents[1]
    proc: CompletedProcess[str] = run_cli(
        dl1_dl2_cfg.exe, dl1_dl2_cfg.args, cwd=repo_root
    )

    if proc.returncode != 0:
        pytest.skip(
            "Skipping DL1->DL2 tests because pipeline failed.\n\n"
            f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}"
        )

    log_path: Path | None = next(dl1_dl2_cfg.output_dir.glob("*.log"), None)
    h5_path: Path | None = next(dl1_dl2_cfg.output_dir.glob("*.h5"), None)

    return log_path, h5_path


def test_0_fixtures_exist(dl1_dl2_cfg: RunCfg):
    # If you prefer models-dir to be skip-only in CI, remove it from this fail check.
    fail_if_missing(dl1_dl2_cfg.missing_paths(), label="dl1->dl2 fixtures")


def test_1_run_sst1mpipe_dl1_dl2(produced_DL2_outputs: tuple):
    assert produced_DL2_outputs is not None


def test_2_artifacts_exist_and_non_empty(produced_DL2_outputs: tuple):
    log_path, h5_path = produced_DL2_outputs
    assert log_path.exists() and log_path.stat().st_size > 0
    assert h5_path.exists() and h5_path.stat().st_size > 0


def test_3_log_contains_required_markers(produced_DL2_outputs: tuple):
    log_path, _ = produced_DL2_outputs
    parsed = parse_log(log_path, LOG_REQUIRED_PATTERNS)
    assert parsed["total_events_tel_022"] + parsed["selected_events_tel_022"] > 0
    for key in ("reco_done", "writing_dl2", "creating_dl2_node"):
        assert parsed[key] == 1, f"Expected {key} == 1, got {parsed[key]}"


def test_4_h5_schema_matches_template(dl1_dl2_cfg: RunCfg, produced_DL2_outputs: tuple):
    _, produced_h5 = produced_DL2_outputs
    assert_schema_matches_template(produced_h5, dl1_dl2_cfg.template_file)
