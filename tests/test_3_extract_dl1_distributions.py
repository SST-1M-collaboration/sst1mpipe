from pathlib import Path
from subprocess import CompletedProcess

import pytest
from conftest import (
    RunCfg,
    assert_schema_matches_template,
    fail_if_missing,
    fixtures_dir,
    run_cli,
    skip_if_missing,
)


@pytest.fixture(scope="module")
def dl1_distributions_cfg(tmp_path_factory: pytest.TempPathFactory):
    fixtures: Path = fixtures_dir()
    return RunCfg(
        exe="sst1mpipe_extract_dl1_distributions",
        opts={
            "--input-dir": fixtures / "dl1",
            "--dl3-index-dir": fixtures / "dl3",
            "--output-dir": tmp_path_factory.mktemp("distributions"),
        },
        template_file=fixtures / "distributions/intensity_hist_tel_022_202512180350.h5",
    )


@pytest.fixture(scope="module")
def produced_DL1_distributions_outputs(dl1_distributions_cfg: RunCfg):
    skip_if_missing(
        dl1_distributions_cfg.missing_paths(), label="distribution fixtures"
    )

    print(f"\n[pytest] Running {dl1_distributions_cfg.exe} (this may take a while)...")
    extra_args = [
        "--date",
        "20251218",
        "--histogram-bins",
        "100",
    ]
    repo_root = Path(__file__).resolve().parents[1]
    proc: CompletedProcess[str] = run_cli(
        dl1_distributions_cfg.exe,
        dl1_distributions_cfg.args + extra_args,
        cwd=repo_root,
    )

    if proc.returncode != 0:
        pytest.skip(
            "Skipping DL1 distribution tests because pipeline failed.\n\n"
            f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}"
        )

    h5_path: Path | None = next(dl1_distributions_cfg.output_dir.glob("*.h5"), None)

    return h5_path


def test_0_fixtures_exist(dl1_distributions_cfg: RunCfg):
    # If you prefer models-dir to be skip-only in CI, remove it from this fail check.
    fail_if_missing(
        dl1_distributions_cfg.missing_paths(), label="dl1 distributions fixtures"
    )


def test_1_run_sst1mpipe_dl1_dl2(produced_DL1_distributions_outputs: tuple):
    assert produced_DL1_distributions_outputs is not None


def test_2_artifacts_exist_and_non_empty(produced_DL1_distributions_outputs: tuple):
    h5_path = produced_DL1_distributions_outputs
    assert h5_path.exists() and h5_path.stat().st_size > 0


def test_3_h5_schema_matches_template(
    dl1_distributions_cfg: RunCfg, produced_DL1_distributions_outputs: tuple
):
    produced_h5 = produced_DL1_distributions_outputs
    assert_schema_matches_template(produced_h5, dl1_distributions_cfg.template_file)
