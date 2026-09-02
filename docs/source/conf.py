# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------

project = "prepmd"
copyright = "2026, CCPBioSim"
author = "CCPBioSim"

version = ""
release = ""


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "nbsphinx",
    "sphinx_copybutton",
    "myst_parser",
]

autosummary_generate = True

autodoc_member_order = "bysource"
autodoc_typehints = "description"

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

templates_path = ["_templates"]
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"

html_theme_options = {
    "light_logo": "logos/prepmd-logo-black-text.svg",
    "dark_logo": "logos/prepmd-logo-white-text.svg",
}

html_static_path = ["_static"]
