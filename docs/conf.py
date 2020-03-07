# -*- coding: utf-8 -*-
import os, sys
from datetime import datetime
#sys.path.insert(0, os.path.abspath('../'))
import biceps

# -- General configuration ------------------------------------------------
extensions = [
    'nbsphinx', # for jupyter notebooks **
    'sphinx.ext.mathjax',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.autosummary',
    'autoapi.sphinx',
]

## Document Python Code
autoapi_type = 'python'
autoapi_dirs = ['biceps'] # Directory of Source Code    (directory of source code)
autoapi_root = 'api'

autodoc_default_flags = ['members', 'inherited-members']
#autodoc_default_flags = ['members']

autosummary_generate = True
autodoc_member_order = 'bysource'
todo_include_todos = False

autoapi_modules = {
        'biceps': {
            'override': False,
            'output': 'auto'
    }
}

numpydoc_class_members_toctree = False

# concatenate both class and __init__ docstrings when generating autodoc class
# docs
autoclass_content = 'both'
#autoclass_content = 'class'

# Execute notebooks before conversion: 'always', 'never', 'auto' (default)
nbsphinx_execute = 'never'

# Use this kernel instead of the one stored in the notebook metadata:
nbsphinx_kernel_name = 'python3'

# If True, the build process is continued even if an exception occurs:
nbsphinx_allow_errors = True

# List of arguments to be passed to the kernel that executes the notebooks:
nbsphinx_execute_arguments = ['--InlineBackend.figure_formats={"png", "pdf"}']

# Controls when a cell will time out (defaults to 30; use -1 for no timeout):
#nbsphinx_timeout = 60

# Default Pygments lexer for syntax highlighting in code cells:
nbsphinx_codecell_lexer = 'ipython3'

# Width of input/output prompts used in CSS:
#nbsphinx_prompt_width = '8ex'

# If window is narrower than this, input/output prompts are on separate lines:
nbsphinx_responsive_width = '250px'

# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = ['.rst', '.md', '.html']

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'BICePs'
authors = u'Yunhui Ge, Robert M. Raddi, Vincent A. Voelz'
date = datetime.now()
copyright = """%s, Temple University, %s\n"""%(date.today().year,authors)
version = biceps.__version__
release = version

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = [
    '_build',
    '_templates',
    '**.ipynb_checkpoints',
    'biceps.egg-info'
]


# This is processed by Jinja2 and inserted before each notebook
nbsphinx_prolog = r"""
{% set docname = env.doc2path(env.docname, base='doc') %}

.. only:: html

    .. role:: raw-html(raw)
        :format: html

    .. nbinfo::

        This page was generated from `{{ docname }}`__.
        Interactive online version:
        :raw-html:`<a href="https://mybinder.org/v2/gh/spatialaudio/nbsphinx/{{ env.config.release }}?filepath={{ docname }}"><img alt="Binder badge" src="https://mybinder.org/badge.svg" style="vertical-align:text-bottom"></a>`

    __ https://github.com/spatialaudio/nbsphinx/blob/
        {{ env.config.release }}/{{ docname }}

.. raw:: latex

    \vfil\penalty-1\vfilneg
    \vspace{\baselineskip}
    \textcolor{gray}{The following section was generated from
    \texttt{\strut{}{{ docname }}}\\[-0.5\baselineskip]
    \noindent\rule{\textwidth}{0.4pt}}
    \vspace{-2\baselineskip}
"""

# This is processed by Jinja2 and inserted after each notebook
nbsphinx_epilog = r"""
.. raw:: latex

    \textcolor{gray}{\noindent\rule{\textwidth}{0.4pt}\\
    \hbox{}\hfill End of
    \texttt{\strut{}{{ env.doc2path(env.docname, base='doc') }}}}
    \vfil\penalty-1\vfilneg
"""

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'BICePsdoc'


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True, #False,
    #'vcs_pageview_mode': '',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Adding our Temp Logo:
html_logo = ''

html_static_path = ['_static']
#html_css_files = [
#        'theme.css'
#]

intersphinx_mapping = {'https://docs.python.org/': None}


# -- Options for Epub output ----------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = authors
epub_publisher = authors
epub_copyright = copyright


#--Adding Markdown-----------------------------------
from recommonmark.parser import CommonMarkParser

source_parsers = {
    '.md': CommonMarkParser,
}

