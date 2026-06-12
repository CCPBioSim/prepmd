# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------

project = "prepmd"
copyright = "2026, CCPBioSim"
author = "CCPBioSim"

version = "1.0.0"
release = "1.0.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "nbsphinx",
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Autosummary / autodoc ---------------------------------------------------

autosummary_generate = True

autosummary_mock_imports = [
    "openff",
    "openff.toolkit",
    "sklearn",
]

autodoc_member_order = "bysource"
autodoc_typehints = "description"

autodoc_mock_imports = [
    "openff",
    "openff.toolkit",
    "sklearn",
]

# -- Napoleon docstring settings --------------------------------------------

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = False
napoleon_use_ivar = True

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_static_path = ["_static"]
