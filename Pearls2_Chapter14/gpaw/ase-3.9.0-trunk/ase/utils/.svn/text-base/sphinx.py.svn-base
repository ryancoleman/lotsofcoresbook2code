from __future__ import print_function
import os
import types
import warnings
from os.path import join
from stat import ST_MTIME
from docutils import nodes
from docutils.parsers.rst.roles import set_classes
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def mol_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    n = []
    t = ''
    while text:
        if text[0] == '_':
            n.append(nodes.Text(t))
            t = ''
            n.append(nodes.subscript(text=text[1]))
            text = text[2:]
        else:
            t += text[0]
            text = text[1:]
    n.append(nodes.Text(t))
    return n, []


def svn_role_tmpl(urlroot,
                  role,
                  rawtext, text, lineno, inliner, options={}, content=[]):
    if text[-1] == '>':
        i = text.index('<')
        name = text[:i - 1]
        text = text[i + 1:-1]
    else:
        name = text
        if name[0] == '~':
            name = name.split('/')[-1]
            text = text[1:]
        if '?' in name:
            name = name[:name.index('?')]
    ref = urlroot + text
    set_classes(options)
    node = nodes.reference(rawtext, name, refuri=ref,
                           **options)
    return [node], []


def trac_role_tmpl(urlroot,
                   role,
                   rawtext, text, lineno, inliner, options={}, content=[]):
    if text[-1] == '>':
        i = text.index('<')
        name = text[:i - 1]
        text = text[i + 1:-1]
    else:
        name = text
        if name[0] == '~':
            name = name.split('/')[-1]
            text = text[1:]
        if '?' in name:
            name = name[:name.index('?')]
    ref = urlroot + text
    set_classes(options)
    node = nodes.reference(rawtext, name, refuri=ref,
                           **options)
    return [node], []


def epydoc_role_tmpl(package_name, urlroot,
                     role,
                     rawtext, text, lineno, inliner, options={}, content=[]):
    name = None
    if text[-1] == '>':
        i = text.index('<')
        name = text[:i - 1]
        text = text[i + 1:-1]

    components = text.split('.')
    if components[0] != package_name:
        components.insert(0, package_name)

    if name is None:
        name = components[-1]

    try:
        module = None
        for n in range(2, len(components) + 1):
            module = __import__('.'.join(components[:n]))
    except ImportError:
        if module is None:
            print('epydoc: could not process:', str(components))
            raise
        for component in components[1:n]:
            module = getattr(module, component)
            ref = '.'.join(components[:n])
            if isinstance(module, (type, types.ClassType)):
                ref += '-class.html'
            else:
                ref += '-module.html'
        if n < len(components):
            ref += '#' + components[-1]
    else:
        ref = '.'.join(components) + '-module.html'

    ref = urlroot + ref
    set_classes(options)
    node = nodes.reference(rawtext, name,
                           refuri=ref,
                           **options)
    return [node], []


def create_png_files(run_all_python_files=False, exclude=[]):
    errcode = os.system('povray -h 2> /dev/null')
    if errcode:
        warnings.warn('No POVRAY!')
        # Replace write_pov with write_png:
        from ase.io import pov
        from ase.io.png import write_png

        def write_pov(filename, atoms, run_povray=False, **parameters):
            p = {}
            for key in ['rotation', 'show_unit_cell', 'radii',
                        'bbox', 'colors', 'scale']:
                if key in parameters:
                    p[key] = parameters[key]
            write_png(filename[:-3] + 'png', atoms, **p)

        pov.write_pov = write_pov

    olddir = os.getcwd()

    for dirpath, dirnames, filenames in os.walk('.'):
        for filename in filenames:
            if filename.endswith('.py'):
                path = join(dirpath, filename)
                if path in exclude:
                    continue
                lines = open(path).readlines()
                try:
                    line = lines[0]
                except IndexError:
                    continue
                if 'coding: utf-8' in line:
                    line = lines[1]
                run = False
                if run_all_python_files:
                    run = True
                elif line.startswith('# creates:'):
                    t0 = os.stat(path)[ST_MTIME]
                    for file in line.split()[2:]:
                        file = file.rstrip(',')
                        try:
                            t = os.stat(join(dirpath, file))[ST_MTIME]
                        except OSError:
                            run = True
                            break
                        else:
                            if t < t0:
                                run = True
                                break
                if run:
                    print('running:', join(dirpath, filename))
                    os.chdir(dirpath)
                    plt.figure()
                    try:
                        execfile(filename, {})
                    finally:
                        os.chdir(olddir)
        if '.svn' in dirnames:
            dirnames.remove('.svn')
