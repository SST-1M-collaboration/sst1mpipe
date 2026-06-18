# Guideline for devellopers
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

## Workflows and Actions 

[ci.yml](.github/workflows/ci.yml)

[documentation.yml](.github/workflows/documentation.yml)

### Building the documentation

The documentation is built using sphynx. All necessary dependencies to build the documentation come with the [environement.yml](environment.yml) file.

```
mamba activate sst1m-dev
cd docs/
sphinx-build . _build/
xdg-open _build/index.html  
```

### Running ruff checks

The code is tested against a set of syntax rules defined in the [project.toml](pyproject.toml) file.
To run these test yourself :

```
mamba activate sst1m-dev
ruff check .
```

### Running the tests

The code contains unit test under the folder [tests/](tests). The notebooks under [notebooks/](notebooks/) are also being tested
To run the tests yourself:

```
mamba activate sst1m-dev
pytest 
```

To see the option used for `pytest`  see the [project.toml](pyproject.toml) file under `[tool.pytest]`
### Running the code test coverage

To see what part of the code in [sst1mpipe/](sst1mpipe) is covered by the tests in `tests/` the `coverage.py`
software can be used as followed :

```
mamba activate sst1m-dev
coverage run
coverage report
```

To see the option used for the code analysis and reporting (i.e. coverage threshold), see the [project.toml](pyproject.toml) file under `[tool.coverage]`