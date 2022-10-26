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

sys.path.insert(0, os.path.abspath("../../"))

import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = "Snekmer"
copyright = (
    "2022 by Pacific Northwest National Laboratory / Battelle Memorial Institute"
)
author = "C. H. Chang, W. C. Nelson, and J. E. McDermott"

# import snekmer

# # The version info for the project you're documenting, acts as replacement for
# # |version| and |release|, also used in various other places throughout the
# # built documents.
# #
# # The short X.Y version.
# version = snekmer.__version__

# if os.environ.get("READTHEDOCS") == "True":
#     # Because Read The Docs modifies conf.py, versioneer gives a "dirty"
#     # version like "5.10.0+0.g28674b1.dirty" that is cleaned here.
#     version = version.partition("+0.g")[0]

# # The full version, including alpha/beta/rc tags.
# release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx_rtd_theme",
    "sphinx_copybutton",
    "sphinxcontrib.bibtex",
    "nbsphinx",
    "sphinxarg.ext",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# autodoc_mock_imports = ["pandas", "snakemake"]

# In lieu of autodoc_mock_imports (not compatible with sphinx-argparse),
# use mock to consolidate pkg requirements for autodoc parsing
import mock

MOCK_MODULES = [
    "Bio",
    "hdbscan",
    "matplotlib",
    "matplotlib.pyplot",
    "numba",
    "numpy",
    "numpy.typing",
    "pandas",
    "scipy",
    "scipy.cluster.hierarchy",
    "scipy.spatial.distance",
    "seaborn",
    "sklearn",
    "sklearn.base",
    "sklearn.cluster",
    "sklearn.decomposition",
    "sklearn.ensemble",
    "sklearn.linear_model",
    "sklearn.manifold",
    "sklearn.metrics",
    "sklearn.metrics.pairwise",
    "sklearn.model_selection",
    "sklearn.pipeline",
    "sklearn.preprocessing",
    "sklearn.svm",
    "sklearn.tree",
    "snakemake",
    "umap",
]
for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = mock.Mock()


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["../_static"]

# set favicon
html_favicon = "../../resources/favicon.svg"

# add logo
html_logo = "../../resources/snekmer_logo.png"
html_theme_options = {
    "logo_only": True,
    "display_version": False,
}

# -- Bibliography settings ---------------------------------------------------
bibtex_bibfiles = ["refs.bib"]
bibtex_encoding = "latin"
bibtex_default_style = "unsrt"

# -- Argparse setup ----------------------------------------------------------
def setup(app):
    app.add_css_file("sphinx-argparse.css")


rst_epilog = f".. |mock_modules| replace:: {MOCK_MODULES}"
