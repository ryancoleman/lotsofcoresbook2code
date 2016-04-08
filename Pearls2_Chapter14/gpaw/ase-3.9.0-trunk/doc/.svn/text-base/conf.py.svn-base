import sys
sys.path.append('.')

extensions = ['ext',
              'images',
              'sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.intersphinx']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'contents'
project = 'ASE'
copyright = '2014, CAMd'

try:
    from ase.version import version
except ImportError:
    version = '3.0.0'
release = version

today_fmt = '%B %d, %Y'
exclude_trees = ['_build']
default_role = 'math'
pygments_style = 'sphinx'
autoclass_content = 'both'
html_style = 'ase.css'
html_logo = '_static/ase.ico'
html_favicon = '_static/ase.ico'
html_static_path = ['_static']
html_last_updated_fmt = '%b %d, %Y'
html_file_suffix = '.html'
htmlhelp_basename = 'ASEdoc'
latex_paper_size = 'a4'
latex_documents = [('contents', 'ase-manual.tex', 'ASE Manual', 'CAMd',
                    'manual')]
latex_preamble = '\usepackage{amsmath}'
intersphinx_mapping = {'gpaw': ('http://wiki.fysik.dtu.dk/gpaw', None),
                       'python': ('http://docs.python.org/2.7', None)}
