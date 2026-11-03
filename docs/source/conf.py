# Configuration file for the Sphinx documentation builder.

import importlib.metadata

# -- Project information -----------------------------------------------------

project = "prepmd"
copyright = "2026, CCPBioSim"
author = "CCPBioSim"

release = importlib.metadata.version("prepmd")


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "myst_parser",
    "autoapi.extension",
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Napoleon docstring settings --------------------------------------------

napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_use_param = True
napoleon_use_ivar = True

# -- API reference generated from the source by sphinx-autoapi --------------

autoapi_dirs = ["../../prepmd"]
autoapi_root = "reference"
autoapi_add_toctree_entry = False
autoapi_member_order = "groupwise"
# Include the __init__ docstring, where the constructor arguments are documented
autoapi_python_class_content = "both"
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
]

# Show type hints in the parameter list rather than the signature
autodoc_typehints = "description"
add_module_names = False
python_maximum_signature_line_length = 88


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "ccpbiosim_sphinx_theme"
html_context = {"github_repo": "prepmd"}
html_static_path = ["_static"]


# Some module variables are assigned twice, such as a default followed by a
# fallback in an except block. Document each name only once.
_documented = set()


def skip_repeated_members(app, what, name, obj, skip, options):
    if name in _documented:
        return True
    _documented.add(name)
    return skip


def setup(app):
    app.connect("autoapi-skip-member", skip_repeated_members)
