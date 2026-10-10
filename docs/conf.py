# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from __future__ import annotations

import os
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from docutils import nodes
    from sphinx import addnodes
    from sphinx import application
    from sphinx import environment

sys.path.insert(0, os.path.abspath('../src'))

from recount3 import __version__

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'recount3'
copyright = '2026, Alexander Alsalihi, Robert M. Flight, Hunter N.B. Moseley'
author = 'Alexander Alsalihi, Robert M. Flight, Hunter N.B. Moseley'

version = __version__
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'

# -- Options for intersphinx extension ---------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#configuration

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
}

# -- Options for todo extension ----------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/todo.html#configuration

todo_include_todos = True

# -- Options for autodoc extension -------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#configuration

autodoc_typehints = 'both'
autoclass_content = 'both'
autodoc_member_order = 'bysource'

_MODULE_ALIASES = {
    'r3_bundle': 'recount3.bundle',
    'r3_config': 'recount3.config',
    'r3_types': 'recount3.types',
}


def _resolve_module_alias(
    app: application.Sphinx,
    env: environment.BuildEnvironment,
    node: addnodes.pending_xref,
    contnode: nodes.TextElement,
) -> nodes.Element | None:
    alias, dot, name = node.get('reftarget', '').partition('.')
    if node.get('refdomain') != 'py' or alias not in _MODULE_ALIASES:
        return None
    if not dot:
        return None
    target = f'{_MODULE_ALIASES[alias]}.{name}'
    if contnode.astext() == node['reftarget']:
        contnode = contnode.__class__(name)
    return env.get_domain('py').resolve_xref(
        env, node['refdoc'], app.builder, node['reftype'], target, node,
        contnode,
    )


def setup(app: application.Sphinx) -> None:
    app.connect('missing-reference', _resolve_module_alias)
