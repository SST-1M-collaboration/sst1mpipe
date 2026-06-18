
# Installation (for developers)
- download latest development version from git repository
- create and activate <b>mamba</b> environment
- install <b>sst1mpipe</b>
```
git clone git@github.com:SST-1M-collaboration/sst1mpipe.git
mamba env create -f sst1mpipe/environment.yml
mamba activate sst1m-dev
pip install -e sst1mpipe
```

# Contributing

# Checking your change 

## Building the documentation

The documentation is built using sphynx. All necessary dependencies to build the documentation come with the [environement.yml](environment.yml) file.

```
mamba activate sst1m-dev
cd docs/
sphinx-build . _build/
xdg-open _build/index.html  
```

## Running ruff checks

The code is tested against a set of syntax rules defined in the [project.toml](pyproject.toml) file.
To run these test yourself :

```
mamba activate sst1m-dev
ruff check .
```

## Running the tests

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