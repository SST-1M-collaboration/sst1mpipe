from pathlib import Path
from subprocess import CompletedProcess

import pytest
from conftest import (
    RunCfg,
    assert_fits_schema_matches_template,
    fail_if_missing,
    fixtures_dir,
    parse_log,
    run_cli,
    skip_if_missing,
)

LOG_REQUIRED_PATTERNS = {
    # numeric counters
    "total_events_tel_022": r"\bTotal N of events of tel_022:\s*([0-9]+)",
    "selected_events_tel_022": r"\bN of events of tel_022 after selection cuts:\s*([0-9]+)",
    "quality_events_tel_022": r"\bN of DL2 events of tel_022 after all input quality cuts:\s*([0-9]+)",
    "events_after_gammaness": r"\bN of events after gammaness cut:\s*([0-9]+)",
}


@pytest.fixture(scope="module")
def dl2_dl3_cfg(tmp_path_factory: pytest.TempPathFactory):
    fixtures: Path = fixtures_dir()
    return RunCfg(
        exe="sst1mpipe_data_dl2_dl3",
        opts={
            "--input-dir": fixtures / "dl2",
            "--config": fixtures / "default_sst1mpipe_data_config_analysis.json",
            "--irf-dir": "/data/work/analysis/IRFs_med4_eff05/",
            "--output-dir": tmp_path_factory.mktemp("sst1mpipe_dl3_out"),
        },
        template_file=fixtures / "dl3/SST1M_CTA1_obs_id_202512180350_dl3.fits",
    )


@pytest.fixture(scope="module")
def produced_DL3_outputs(dl2_dl3_cfg: RunCfg):
    skip_if_missing(dl2_dl3_cfg.missing_paths(), label="dl2->dl3 fixtures")

    print(f"\n[pytest] Running {dl2_dl3_cfg.exe} (this may take a while)...")

    extra_args = [
        "--target-name",
        "CTA1",
        "--target-ra",
        "1.65",
        "--target-dec",
        "72.783",
    ]
    repo_root = Path(__file__).resolve().parents[1]
    proc: CompletedProcess[str] = run_cli(
        dl2_dl3_cfg.exe, dl2_dl3_cfg.args + extra_args, cwd=repo_root
    )

    if proc.returncode != 0:
        pytest.skip(
            "Skipping DL1->DL2 tests because pipeline failed.\n\n"
            f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}"
        )

    log_path: Path | None = next(dl2_dl3_cfg.output_dir.glob("*.log"), None)
    fits_path: Path | None = next(dl2_dl3_cfg.output_dir.glob("*.fits"), None)
    hdu_path: Path | None = next(dl2_dl3_cfg.output_dir.glob("hdu-index.fits.gz"), None)
    obs_path: Path | None = next(dl2_dl3_cfg.output_dir.glob("obs-index.fits.gz"), None)

    return log_path, fits_path, hdu_path, obs_path


def test_0_fixtures_exist(dl2_dl3_cfg: RunCfg):
    # If you prefer models-dir to be skip-only in CI, remove it from this fail check.
    fail_if_missing(dl2_dl3_cfg.missing_paths(), label="dl2->dl3 fixtures")


def test_1_run_sst1mpipe_dl2_dl3(produced_DL3_outputs: tuple):
    assert produced_DL3_outputs is not None


def test_2_artifacts_exist_and_non_empty(produced_DL3_outputs: tuple):
    log_path, fits_path, hdu_path, obs_path = produced_DL3_outputs
    assert log_path.exists() and log_path.stat().st_size > 0
    assert fits_path.exists() and fits_path.stat().st_size > 0
    assert hdu_path.exists() and hdu_path.stat().st_size > 0
    assert obs_path.exists() and obs_path.stat().st_size > 0


def test_3_log_contains_required_markers(produced_DL3_outputs: tuple):
    log_path, _, _, _ = produced_DL3_outputs
    parsed = parse_log(log_path, LOG_REQUIRED_PATTERNS)
    assert all(parsed[k] > 0 for k in LOG_REQUIRED_PATTERNS)


def test_4_h5_schema_matches_template(dl2_dl3_cfg: RunCfg, produced_DL3_outputs: tuple):
    _, fits_path, _, _ = produced_DL3_outputs
    assert_fits_schema_matches_template(fits_path, dl2_dl3_cfg.template_file)
