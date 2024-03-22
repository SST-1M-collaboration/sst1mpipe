# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

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
    "sphinx_automodapi.automodapi",
]

numpydoc_show_class_members = False
autosummary_generate = True

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3.8", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "astropy": ("https://docs.astropy.org/en/latest/", None),
    "pytables": ("https://www.pytables.org/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "matplotlib": ("https://matplotlib.org/", None),
    "traitlets": ("https://traitlets.readthedocs.io/en/stable/", None),
    "ctapipe": ("https://ctapipe.readthedocs.io/en/v0.17.1/", None)
}

# These links are ignored in the checks, necessary due to broken intersphinx for
# these
nitpick_ignore = [
    ("py:class", "sst1mpipe.instrument.camera.Camera"),
    ("py:class", "cts_core.camera.Camera"),
    ("py:class", "ctapipe.instrument.camera.geometry.CameraGeometry"),
    ("py:class", "ctapipe.core.tool.Tool"),
    ("py:class", "ctapipe.core.component.Component"),
    ("py:class", "ctapipe.core.container.Container"),
    ("py:class", "ctapipe.calib.camera.flatfield.FlatFieldCalculator"),
    ("py:class", "ctapipe.calib.camera.pedestals.PedestalCalculator"),
    # coming from inherited traitlets docs
    ("py:class", "t.Union"),
    ("py:class", "t.Dict"),
    ("py:class", "t.Tuple"),
    ("py:class", "t.List"),
    ("py:class", "t.Any"),
    ("py:class", "t.Type"),
    ("py:class", "Config"),
    ("py:class", "Unicode"),
    ("py:class", "StrDict"),
    ("py:class", "ClassesType")
]


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
import sphinx_rtd_theme
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']