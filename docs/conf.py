# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import shutil
from pathlib import Path

import sst1mpipe

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'sst1mpipe'
copyright = '2024, SST-1M collaboration'
author = 'SST-1M collaboration'
release = sst1mpipe.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    'sphinx.ext.githubpages',
    "numpydoc",
    "sphinx_rtd_theme",
    "sphinx_togglebutton",
    "nbsphinx",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Tutorials (notebooks) ---------------------------------------------------
# The notebooks live in notebooks/ (outside of the Sphinx source directory, and run by the
# tests with nbmake): the ones of the documentation are copied to docs/tutorials/ at build time.
# They are rendered with their saved outputs, not executed (the docs build has no test data).
TUTORIAL_NOTEBOOKS = [
    "sst1m_read_zfits_event_source.ipynb",
]
_notebooks_dir = Path(__file__).parent.parent / "notebooks"
_tutorials_dir = Path(__file__).parent / "tutorials"
_tutorials_dir.mkdir(exist_ok=True)
for _notebook in TUTORIAL_NOTEBOOKS:
    shutil.copy2(_notebooks_dir / _notebook, _tutorials_dir / _notebook)

nbsphinx_execute = "never"
nbsphinx_prolog = """
.. note::

   This page is the notebook ``notebooks/{{ env.docname.split('/')[-1] }}.ipynb`` of the
   `sst1mpipe repository <https://github.com/SST-1M-collaboration/sst1mpipe>`_:
   download it to run it yourself.
"""



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
