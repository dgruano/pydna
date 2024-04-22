# -*- coding: utf-8 -*-
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pydna"
copyright = "2024, Björn F. Johansson"
author = "Björn F. Johansson"
import os
import sys
sys.path.insert(0, os.path.abspath("../src/pydna"))
# contents of docs/conf.py
from importlib.metadata import version

release = version("pydna")
# for example take major/minor
version = ".".join(release.split(".")[:3])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx_rtd_theme",
]

# Add mappings https://kev.inburke.com/kevin/sphinx-interlinks
intersphinx_mapping = {
    "biopython": ("https://biopython.org/docs/latest/api/", None),
    "python": ("http://docs.python.org/3.8", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_member_order = 'bysource'
autodoc_preserve_defaults = True

numpydoc_show_class_members = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

texinfo_documents = [
    (
        "index",
        "Pydna",
        "Pydna Documentation",
        "Björn Johansson",
        "Pydna",
        "One line description of project.",
        "Miscellaneous",
    ),
]
