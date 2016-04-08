#!/usr/bin/python
import os
import sys
import glob
import tempfile

if '--dir' in sys.argv:
    i = sys.argv.index('--dir')
    dir = sys.argv[i + 1]
else:
    dir = None

if '--email' in sys.argv:
    i = sys.argv.index('--email')
    email = sys.argv[i + 1]
else:
    email = None

if '--tarfiledir' in sys.argv:
    i = sys.argv.index('--tarfiledir')
    tarfiledir = sys.argv[i + 1]
else:
    tarfiledir = None

tmpdir = tempfile.mkdtemp(prefix='gpaw-sphinx-', dir=dir)
os.chdir(tmpdir)


def build():
    if os.system('svn export ' +
                 'https://svn.fysik.dtu.dk/projects/ase/trunk ase') != 0:
        raise RuntimeError('Checkout of ASE failed!')
    os.chdir('ase')
    if os.system('python setup.py install --home=..') != 0:
        raise RuntimeError('Installation of ASE failed!')
    os.chdir('..')
    if os.system('svn checkout ' +
                 'https://svn.fysik.dtu.dk/projects/gpaw/trunk gpaw') != 0:
        raise RuntimeError('Checkout of GPAW failed!')
    os.chdir('gpaw')
    if os.system('python setup.py install --home=.. 2> error') != 0:
        raise RuntimeError('Installation of GPAW failed!')

    assert os.system('wget --no-check-certificate --quiet ' +
                     'http://wiki.fysik.dtu.dk/gpaw-files/' +
                     'gpaw-setups-latest.tar.gz') == 0

    assert os.system('tar xvzf gpaw-setups-latest.tar.gz') == 0

    setups = tmpdir + '/gpaw/' + glob.glob('gpaw-setups-[0-9]*')[0]

    # Generate tar-file:
    assert os.system('python setup.py sdist') == 0

    if os.system('epydoc --docformat restructuredtext --parse-only ' +
                 '--name GPAW ' +
                 '--url http://wiki.fysik.dtu.dk/gpaw ' +
                 '--show-imports --no-frames -v gpaw &> epydoc.out') != 0:
        raise RuntimeError('Epydoc failed!')

    epydoc_output = open('epydoc.out').readlines()
    errors = []
    for line in epydoc_output:
        if line[0] == '|' or line[:2] == '+-':
            errors.append(line)
    if errors:
        fd = open('epydoc.errors', 'w')
        fd.write(''.join(errors))
        fd.close()
        if email is not None:
            assert os.system(
                'mail -s "GPAW: EpyDoc errors" %s < epydoc.errors' %
                email) == 0

    # ase installs under lib independently of the platform
    sys.path.insert(0, '%s/lib/python' % tmpdir)
    # gpaw installs under libdir
    from distutils.sysconfig import get_config_var
    libdir = os.path.split(get_config_var('LIBDIR'))[-1]
    # import gpaw from where it was installed
    sys.path.insert(0, '%s/%s/python' % (tmpdir, libdir))

    from gpaw.version import version

    os.chdir('doc')
    assert os.system('sed -i s/gpaw-snapshot/gpaw-%s/ download.rst' %
                     version) == 0
    os.mkdir('_build')
    if os.system('PYTHONPATH=%s/%s/python:%s/lib/python ' %
                 (tmpdir, libdir, tmpdir) +
                 'GPAW_SETUP_PATH=%s ' % setups +
                 'sphinx-build -n . _build') != 0:
        raise RuntimeError('Sphinx failed!')
    assert os.system('cd _build; cp _static/searchtools.js .') == 0
    # Don't serve mixed content:
    assert os.system(
        'cd _build; find -name "*.html" | xargs sed -i ' +
        '"s%http://cdn.mathjax.org%https://cdn.mathjax.org%"') == 0
    if 0:
        if os.system('sphinx-build -b latex . _build') != 0:
            raise RuntimeError('Sphinx failed!')
        os.chdir('_build')
        if os.system('make gpaw-manual.pdf') != 0:
            raise RuntimeError('pdflatex failed!')
    else:
        os.chdir('_build')

    assert os.system('mv ../../html epydoc;' +
                     'mv ../../dist/gpaw-%s.tar.gz .' % version) == 0

if tarfiledir is not None:
    try:
        os.remove(tarfiledir + '/gpaw-webpages.tar.gz')
    except OSError:
        pass

build()
    
if tarfiledir is not None:
    os.system('cd ..; tar czf %s/gpaw-webpages.tar.gz _build' % tarfiledir)

os.system('cd; /bin/rm -rf ' + tmpdir)
