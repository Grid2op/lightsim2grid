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
release = "0.11.0.dev0"
version = '0.10'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.imgmath',
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

    # markdown renderer
    # 'm2r',
    'recommonmark',
]

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
html_experimental_html5_writer = True
html_theme = "sphinx_rtd_theme" #"alabaster" #'basic' # 'alabaster'
highlight_language = 'python3'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

