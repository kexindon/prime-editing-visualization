#Sphinx build for the visualizer's user guide.
#
#The page content lives in DOCUMENTATION.md at the repository root, so that one
#file serves both GitHub (where people land first) and this site. MyST parses
#the Markdown; a symlink in this directory points at it.

project = 'pegRNA Design Visualizer'
copyright = '2026'
author = 'Kexin Dong'

extensions = ['myst_parser']

#Match the PEGG 3 documentation, which this tool is a companion to.
html_theme = 'nltk_theme'

source_suffix = {'.md': 'markdown', '.rst': 'restructuredtext'}
master_doc = 'index'

#The images sit in docs/media/, referenced from DOCUMENTATION.md by paths
#relative to the repository root. Copying that tree in keeps one set of paths
#working for both GitHub and this build.
html_static_path = []
exclude_patterns = ['_build', 'docs/media/README.md']

myst_enable_extensions = ['linkify', 'colon_fence']
myst_heading_anchors = 3
