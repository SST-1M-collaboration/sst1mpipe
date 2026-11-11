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

# `name` in the docstrings is a reference to the Python object
default_role = "py:obj"

# members documented with autodoc :members:, not in tables (they would list all the inherited
# traitlets methods of the ctapipe components)
numpydoc_show_class_members = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "ctapipe": ("https://ctapipe.readthedocs.io/en/stable/", None),
}
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Tutorials (notebooks) ---------------------------------------------------
# The notebooks live in notebooks/ (outside of the Sphinx source directory, and run by the
# tests with nbmake): the ones of the documentation are copied to docs/tutorials/ at build time.
# They are rendered with their saved outputs, not executed (the docs build has no test data).
TUTORIAL_NOTEBOOKS = [
    "sst1m_read_zfits_event_source.ipynb",
    "sst1m_trigger_traces.ipynb",
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


# -- Configuration options of the ctapipe components -------------------------
# The configuration options (traits with config=True) of the documented components are listed as
# attributes of the class, with their type, default and help; the traits are not documented again
# as members.
def _is_trait(obj):
    from traitlets import TraitType

    return isinstance(obj, TraitType)


def _add_configuration_options(app, what, name, obj, options, lines):
    from traitlets.config import Configurable

    if what != "class" or not (isinstance(obj, type) and issubclass(obj, Configurable)):
        return
    traits = obj.class_traits(config=True)
    if not traits:
        return
    lines += ["", ".. rubric:: Configuration options", ""]
    for trait_name, trait in sorted(traits.items()):
        help_text = " ".join(trait.help.split()) or "(no description)"
        lines += [
            f".. py:attribute:: {trait_name}",
            f"   :type: {type(trait).__name__}",
            f"   :value: {trait.default_value_repr()}",
            "",
            f"   {help_text}",
            "",
        ]


def _skip_traits(app, what, name, obj, skip, options):
    return True if _is_trait(obj) else None


def setup(app):
    # after numpydoc, which reformats the docstrings
    app.connect("autodoc-process-docstring", _add_configuration_options, priority=600)
    app.connect("autodoc-skip-member", _skip_traits)
