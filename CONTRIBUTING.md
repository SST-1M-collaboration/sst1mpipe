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

