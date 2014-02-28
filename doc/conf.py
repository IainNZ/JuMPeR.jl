# -*- coding: utf-8 -*-
import sys, os, juliadoc

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('sphinx'))

# -- General configuration -----------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx.ext.mathjax', 'juliadoc.julia', 'juliadoc.jlhelp']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'jumper'

# General information about the project.
project = u'JuMPeR -- Julia for Mathematical Programming'
AUTHORS = u"Iain Dunning"
copyright = u'2014, '+AUTHORS

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '0.0'
# The full version, including alpha/beta/rc tags.
release = '0.0'

exclude_patterns = ['_build']
pygments_style = 'sphinx'
primary_domain = 'jl'
highlight_language = 'julia'


# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'julia'
html_theme_path = [juliadoc.get_theme_dir()]
html_sidebars = juliadoc.default_sidebars()
htmlhelp_basename = 'JuMPeRJlDoc'


# -- Options for LaTeX output --------------------------------------------------

latex_elements = {
    'utf8extra': r'''
        \DeclareUnicodeCharacter{00A0}{\nobreakspace}
        \DeclareUnicodeCharacter{2203}{\ensuremath{\exists}}
        \DeclareUnicodeCharacter{2200}{\ensuremath{\forall}}
        \DeclareUnicodeCharacter{27FA}{\ensuremath{\Longleftrightarrow}}
    ''',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
  ('jumper', 'JuMPeR.jl.tex', u'JuMPeR.jl Documentation',
   AUTHORS, 'manual'),
]


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('jumper', 'JuMPeR', u'JuMPeR.jl Documentation',
     [AUTHORS], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output ------------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('jumper', 'JuMPeR', u'JuMPeR.jl Documentation',
   AUTHORS, 'JuMPeR', 'JUMP - Julia for Mathematical Programming.',
   'Miscellaneous'),
]