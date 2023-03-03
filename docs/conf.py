# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
<<<<<<< HEAD
# https://www.sphinx-doc.org/en/master/usage/configuration.html
=======
# http://www.sphinx-doc.org/en/master/config
>>>>>>> chlamdb/master

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

<<<<<<< HEAD
project = 'annotation'
copyright = '2020, Trestan Pillonel, Bastian Marquis'
author = 'Trestan Pillonel, Bastian Marquis'

# The full version, including alpha/beta/rc tags
release = '1.0'

=======
project = 'ChlamDB'
copyright = '2019, Trestan Pillonel'
author = 'Trestan Pillonel'

# The full version, including alpha/beta/rc tags
release = '2.0 (June 2019)'

todo_include_todos = True
>>>>>>> chlamdb/master

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
<<<<<<< HEAD
=======
    'sphinx.ext.todo',
>>>>>>> chlamdb/master
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
<<<<<<< HEAD
exclude_patterns = []
=======
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
>>>>>>> chlamdb/master


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
<<<<<<< HEAD
html_theme = 'nature'
=======

import quark_sphinx_theme
import sphinx_bootstrap_theme
# import sphinx_theme_material
html_theme_path = [quark_sphinx_theme.get_path()]
html_theme = 'nature'
#html_theme = 'bootstrap'
#html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()
#html_theme = 'alabaster'
>>>>>>> chlamdb/master

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
<<<<<<< HEAD
=======
html_css_files = [
    'css/bootstrap.min.css',
    'css/custom.css',
]

html_js_files = [
    'js/bootstrap.min.js',
    'js/custom.js'
]

>>>>>>> chlamdb/master
