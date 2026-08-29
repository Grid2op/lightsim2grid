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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'LightSim2Grid'
copyright = '2019, RTE France'
author = 'Benjamin DONNOT'

# The full version, including alpha/beta/rc tags
release = "1.0.0"
version = '1.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',  # client-side math rendering: no local LaTeX install needed
                           # (unlike sphinx.ext.imgmath, which shells out to `latex`)
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    # 'builder',
    'sphinx.ext.extlinks',
    'sphinx.ext.napoleon',
    'sphinxcontrib_trio',
     "sphinx_rtd_theme",
    # toc of modules
    'autodocsumm',
    # 'sphinx.ext.autosectionlabel',

    # 'details',
    #'exception_hierarchy',

    # for pdf
    # 'rst2pdf.pdfbuilder'
]

# External inventories so :class:/:func:/:attr: references to third-party types
# used throughout the docstrings (numpy/scipy/pandas type hints, grid2op and
# pypowsybl classes) resolve instead of showing up as nitpicky "reference target
# not found" warnings. Fetched over the network at build time (and cached) --
# Read the Docs builds have normal internet access; a fully offline build just
# falls back to unresolved (but non-fatal) references for these.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'grid2op': ('https://grid2op.readthedocs.io/en/latest/', None),
    'pypowsybl': ('https://pypowsybl.readthedocs.io/en/latest/', None),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = []  # ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme" #"alabaster" #'basic' # 'alabaster'
highlight_language = 'python3'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

