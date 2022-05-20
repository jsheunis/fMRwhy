# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

# sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath("sphinxext"))
sys.path.insert(0, os.path.abspath(os.path.pardir))
this_dir = os.path.dirname(os.path.abspath(__file__))
matlab_src_dir = os.path.abspath(os.path.join(this_dir, ".."))
primary_domain = "mat"


# -- Project information -----------------------------------------------------

project = "fMRwhy"
copyright = "2020, J.S. Heunis"
author = "J.S. Heunis"

# The full version, including alpha/beta/rc tags
release = "v0.0.2"


# -- General configuration ---------------------------------------------------


# generate autosummary even if no references
# autosummary_generate = True

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.viewcode",
    "sphinxcontrib.matlab",
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# installing theme package
html_theme = "sphinx_rtd_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.

html_theme_options = {
    "includehidden": False,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
