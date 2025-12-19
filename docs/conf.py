import os

# import subprocess
import sys

sys.path.insert(0, os.path.abspath(".."))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "fulldia"
copyright = "2025, sj"
author = "sj"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
]
napoleon_numpy_docstring = True
napoleon_google_docstring = False

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = ["_static"]

# api_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "api_reference"))
# if not os.path.exists(api_dir):
#     os.makedirs(api_dir)
#
# subprocess.run(
#     [
#         "sphinx-apidoc",
#         "-o",
#         api_dir,
#         os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "full_dia")),
#         "--force",
#         "--separate",
#     ],
#     check=True,
# )
