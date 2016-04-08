from ase.utils.sphinx import mol_role
from ase.utils.sphinx import svn_role_tmpl, trac_role_tmpl, epydoc_role_tmpl
from ase.utils.sphinx import create_png_files


def svn_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    return svn_role_tmpl('http://svn.fysik.dtu.dk/projects/'
                         'ase/trunk/',
                         role,
                         rawtext, text, lineno, inliner, options, content)


def trac_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    return trac_role_tmpl('http://trac.fysik.dtu.dk/projects/'
                          'ase/browser/trunk/',
                          role,
                          rawtext, text, lineno, inliner, options, content)


def epydoc_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    return epydoc_role_tmpl('ase', 'http://wiki.fysik.dtu.dk/ase/epydoc/',
                            role,
                            rawtext, text, lineno, inliner, options, content)


def setup(app):
    app.add_role('mol', mol_role)
    app.add_role('svn', svn_role)
    app.add_role('trac', trac_role)
    app.add_role('epydoc', epydoc_role)
    #import atexit
    #atexit.register(fix_sidebar)
    create_png_files()
