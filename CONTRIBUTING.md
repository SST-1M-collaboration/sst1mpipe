# Guideline for developers

## Installation

To set up a development environment, follow these steps:

1. **Clone the repository:**

```bash
git clone git@github.com:SST-1M-collaboration/sst1mpipe.git
```

2. **Create the mamba environment** using the provided [environment file](environment.yml):

```bash
mamba env create -f sst1mpipe/environment.yml
```

3. **Activate the environment:**

```bash
mamba activate sst1m-dev
```

4. **Install `sst1mpipe` in editable mode:**

```bash
pip install -e sst1mpipe
```

Installing the package in editable mode (`-e`) allows your local modifications to the source code to be immediately reflected without reinstalling the package after each change.

## Contributing

Thank you for your interest in contributing to this project!

The typical workflow for submitting changes is the following:

1. **Open an issue** describing the feature you would like to implement or the bug you would like to fix.

2. **Create a new branch** from the `main` branch:

```bash
git checkout -b branch-name
```

3. **Implement your changes** and commit them incrementally.

It is recommended to create **small, focused commits** that each correspond to a logical change (for example, fixing a bug, adding a feature, or updating documentation). This makes it easier to review your work, understand the history of the project, and revert specific changes if needed.

```bash
git commit -m "A short message describing the change"
```

You are encouraged to commit your work regularly as you progress rather than creating a single large commit at the end.


4. **Push your branch** to the GitHub repository (`origin`):

```bash
git push origin branch-name
```

5. **Open a Pull Request (PR)** and provide a clear description of the changes you have implemented.

6. **Ensure that all tests pass** and fix any issues if necessary.

The project uses automated GitHub Actions workflows to run tests, check code quality, and verify that changes do not introduce regressions. After opening your Pull Request, make sure that all checks complete successfully and address any reported issues if necessary.

For more information about the available workflows and how to run them locally, see the [Workflows and Actions](#workflows-and-actions) section below.

7. **Request a review** and address any feedback or requested changes from the reviewers.

8. 🎉 **Celebrate your contribution to improving the project!**

## GitHub Workflows and Actions

GitHub Workflows (via **GitHub Actions**) are powerful automation tools available in GitHub repositories. They are used to automatically run jobs such as testing, linting, and documentation builds whenever changes are pushed or pull requests are opened.

All workflows for this project are defined in the [.github/worflows/](.github/workflows/)
folder.

Key workflows include:

* [ci.yml](.github/workflows/ci.yml) — Continuous integration (tests, linting, and quality checks)
* [documentation.yml](.github/workflows/documentation.yml) — Builds and deploys the project documentation
* [release.yml](.github/workflows/release.yml) — Builds the package and publishes it on PyPI when a version tag is pushed, see [Making a release](#making-a-release)

---

## Project Files Overview

This project relies on a few key files that define how the codebase is tested, formatted, and managed:

### `pyproject.toml`

The [pyproject.toml](pyproject.toml) file is the central configuration file for most Python tooling in this project.

It is used to configure:

* **Setuptools** (how to install the software, create command line applications)
* **Ruff** (linting and code style rules)
* **Pytest** (test configuration and options)
* **Coverage.py** (test coverage thresholds and reporting)
* Potentially other tools used in the development workflow

### `MANIFEST.in`

The [MANIFEST.in](MANIFEST.in) file selects the files of the source distribution (`.tar.gz`) published on
PyPI. It removes the test data (`sst1mpipe/resources/test_data`, ~550 MB of Git LFS files) from it; the
wheel excludes them in `pyproject.toml`.

### `environment.yml`

### `.gitattributes`

The `.gitattributes` file defines how Git handles specific file types in the repository.
It is especially important for:

* Ensuring consistent line endings across platforms
* Configuring Git LFS tracking for large files
* Defining diff/merge behavior for specific file types

### `.gitignore`

### `.readthedocs.yaml`

### `LICENSE`

### `README.md`

### `CONTRIBUTING.md`

## Building the Documentation

The project documentation is built using **Sphinx**.
To build and preview the documentation locally:

```bash
sphinx-build docs/ _build/
xdg-open _build/index.html
```
## Running Ruff (Code Style Checks)

The codebase is checked against a set of style and linting rules defined in [pyproject.toml](pyproject.toml).

To run the checks locally:

```bash
ruff check 
```

## Running the Tests

Unit tests are located in the [tests/](tests/) directory. In addition, Jupyter notebooks in [notebooks/](notebooks/) are included in the test suite.
Configuration options for `pytest` can be found in [pyproject.toml](pyproject.toml) under the `[tool.pytest]` section.

To run the full test suite:

```bash
pytest 
```

## Code Coverage

To measure which parts of the codebase are covered by tests, the project uses `coverage.py`.
Coverage thresholds and reporting options are defined in [pyproject.toml](pyproject.toml) under the `[tool.coverage]` section.

Run coverage analysis and show the result with:

```bash
coverage report
```

## Working with Large Files

This repository uses **Git LFS (Large File Storage)** to handle large files that should not be stored directly in Git history.

Make sure Git LFS is installed and initialized:

```bash
git lfs install
```

Large file tracking rules are defined in `.gitattributes`, which specifies which file types should be handled by Git LFS instead of standard Git storage.

## Making a release

How to publish a new version of `sst1mpipe` on [PyPI](https://pypi.org/project/sst1mpipe/).

### Version numbers

Since 0.32, the versions of `sst1mpipe` follow the versions of
[ctapipe](https://github.com/cta-observatory/ctapipe): `sst1mpipe` 0.32.x uses ctapipe 0.32.
The versions follow [PEP 440](https://peps.python.org/pep-0440/):

| Version     | Meaning                                                                 |
|-------------|-------------------------------------------------------------------------|
| `0.32.0a1`  | alpha release (`a2`, `a3`, … for the next alpha releases)              |
| `0.32.0b1`  | beta release                                                            |
| `0.32.0rc1` | release candidate                                                       |
| `0.32.0`    | final release                                                           |
| `0.32.1`    | bug fix release, still with ctapipe 0.32                                |
| `0.33.0`    | first release with ctapipe 0.33 (update the ctapipe pin in `pyproject.toml` and `environment.yml`) |

The 0.32.0aN versions are alpha releases: the R0 to DL1 processing (`sst1mpipe-monitoring`,
`sst1mpipe-process`) uses the new ctapipe tools, while the DL1 to DL3 steps and the stereo matching are
not migrated yet.

`pip install sst1mpipe` ignores the alpha, beta and rc versions: they are installed with
`pip install --pre sst1mpipe` or with the exact version, e.g. `pip install sst1mpipe==0.32.0a1`.

> [!WARNING]
> A version published on PyPI can never be uploaded again, even after deleting it: check the
> package before pushing the tag (see [Checking the package](#2-check-the-package-locally)).

### How a release works

The workflow [.github/workflows/release.yml](.github/workflows/release.yml) runs when a tag `v*` is
pushed:

1. **build**: checks that the tag is `v` + the version of `pyproject.toml` (tag `v0.32.0a2` for
   version `0.32.0a2`), builds the source distribution and the wheel (`python -m build`) and checks
   them (`twine check --strict`).
2. **publish-pypi**: uploads them on PyPI with
   [trusted publishing](https://docs.pypi.org/trusted-publishers/): PyPI accepts the upload because it
   comes from this workflow of this repository, without any token or password.

The workflow can also be started by hand (Actions → release → Run workflow): it then only builds and
checks the package, without publishing it.

The test data (`sst1mpipe/resources/test_data`, ~550 MB of Git LFS files) are not in the package (PyPI
limit: 100 MB per file): they are excluded in `pyproject.toml` (wheel) and `MANIFEST.in` (source
distribution). The tests are run from the git repository.

### Release steps

#### 1. Prepare the version

In a branch and a pull request:

- set the new `version` in [pyproject.toml](pyproject.toml);
- if the version follows a new ctapipe version, update the ctapipe pin in `pyproject.toml`
  (`ctapipe[eventio]~=0.33.0`) and in `environment.yml`;
- check that the dependencies of `pyproject.toml` match `environment.yml`;
- update the installation instructions of [README.md](README.md) and `docs/introduction.rst` if needed.

#### 2. Check the package locally

```bash
pip install build twine
python -m build
twine check --strict dist/*
```

Check the size and the content of the wheel (no test data, ~2 MB), then install it in a new
environment and try it:

```bash
unzip -l dist/sst1mpipe-*.whl
python -m venv /tmp/sst1mpipe-release
/tmp/sst1mpipe-release/bin/pip install dist/sst1mpipe-*.whl
/tmp/sst1mpipe-release/bin/sst1mpipe-process --help
```

#### 3. Merge and tag

Once the pull request is merged and the CI of `main` passes, tag the merge commit and push the tag:

```bash
git checkout main
git pull origin main
git tag v0.32.0a2
git push origin v0.32.0a2
```

Follow the workflow in the Actions tab of GitHub. Once it is done, the version is on
https://pypi.org/project/sst1mpipe/.

#### 4. Check the release

```bash
python -m venv /tmp/sst1mpipe-check
/tmp/sst1mpipe-check/bin/pip install --pre sst1mpipe==0.32.0a2
/tmp/sst1mpipe-check/bin/python -c "import sst1mpipe; print(sst1mpipe.__version__)"
```

#### 5. Publish the GitHub release

On GitHub, Releases → Draft a new release, choose the tag, generate the release notes (list of the
merged pull requests) and add a short summary of the main changes. Tick "Set as a pre-release" for
the alpha, beta and rc versions.

### Troubleshooting

| Error in the workflow | Cause and fix |
|-----------------------|---------------|
| `tag v… is not the version v… of pyproject.toml` | The tag does not match the version. Nothing is published: delete the tag (`git tag -d vX` and `git push origin :refs/tags/vX`), fix the version or the tag, and push the tag again. |
| `invalid-publisher: valid token, but no corresponding publisher` | The trusted publisher is missing or different on PyPI, see [One-time setup](#one-time-setup-of-pypi). |
| `File already exists` / `400 Bad Request` | This version is already on PyPI and cannot be uploaded again: release a new version. |
| `Failed to parse pyproject.toml` | Invalid TOML in `pyproject.toml` (e.g. a missing comma): check it locally with `python -m build`. |

If a job fails after the build (e.g. the upload), fix the cause and use "Re-run failed jobs" in the
Actions tab: the same tag can be used again as long as nothing was uploaded.

### One-time setup of PyPI

Needed once per package (`sst1mpipe`, `sst1m-cts`), before its first release.

1. Create an account on https://pypi.org and enable two-factor authentication (required to publish).
2. In Account settings → Publishing, add a (pending) trusted publisher of type GitHub:

   | Field             | Value                  |
   |-------------------|------------------------|
   | PyPI Project Name | `sst1mpipe`            |
   | Owner             | `SST-1M-collaboration` |
   | Repository name   | `sst1mpipe`            |
   | Workflow name     | `release.yml`          |
   | Environment name  | `pypi`                 |

3. After the first release, add the other maintainers in the project on PyPI (Manage → Collaborators),
   so that the package does not depend on a single account.

### Dependency: sst1m-cts

`sst1mpipe` uses the camera configuration of DigiCam of [CTS](https://github.com/SST-1M-collaboration/CTS),
published on PyPI as [sst1m-cts](https://pypi.org/project/sst1m-cts/) (import `sst1m_cts_core`).
CTS is released the same way (tag `v*`, workflow `release.yml` of the CTS repository). A new version of
`sst1mpipe` that needs a new version of CTS requires the CTS release first, then the update of
`sst1m-cts>=…` in `pyproject.toml` and `environment.yml`.
