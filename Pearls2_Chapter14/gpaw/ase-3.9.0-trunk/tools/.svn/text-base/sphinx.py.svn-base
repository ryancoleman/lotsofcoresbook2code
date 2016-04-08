#!/usr/bin/python
import os
import sys
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

if '--nodoc' in sys.argv:
    i = sys.argv.index('--nodoc')
    doc = False
else:
    doc = True

if '--tarfiledir' in sys.argv:
    i = sys.argv.index('--tarfiledir')
    tarfiledir = sys.argv[i + 1]
else:
    tarfiledir = None

tmpdir = tempfile.mkdtemp(prefix='ase-sphinx-', dir=dir)
os.chdir(tmpdir)


def build(email):
    if os.system('svn checkout ' +
                 'https://svn.fysik.dtu.dk/projects/ase/trunk ase') != 0:
        raise RuntimeError('Checkout of ASE failed!')
    os.chdir('ase')
    if os.system('python setup.py install --home=..') != 0:
        raise RuntimeError('Installation failed!')

    # ase installs under lib independently of the platform
    sys.path.insert(0, '%s/lib/python' % tmpdir)

    from ase.test import test
    from ase.version import version

    if 'ASE_CALCULATORS' in os.environ:
        calculators = os.environ['ASE_CALCULATORS'].split(',')
    else:
        calculators = []

    # Run test-suite:
    stream = open('test-results.txt', 'w')
    results = test(verbosity=2, calculators=calculators,
                   display=False, stream=stream)
    stream.close()
    if len(results.failures) > 0 or len(results.errors) > 0:
        if email is not None:
            subject = 'ASE %s: test-suite failed!' % str(version)
            os.system('mail -s "%s" %s < %s' %
                      (subject, email, 'test-results.txt'))
        raise RuntimeError('Testsuite failed!')

    # Generate tar-file:
    assert os.system('python setup.py sdist') == 0

    if doc:

        if os.system('epydoc --docformat restructuredtext --parse-only ' +
                     '--name ASE ' +
                     '--url http://wiki.fysik.dtu.dk/ase ' +
                     '--show-imports --no-frames -v ase &> epydoc.out') != 0:
            raise RuntimeError('Epydoc failed!')

        epydoc_errors = open('epydoc.out').read()
        if ' Warning:' in epydoc_errors:
            sys.stderr.write(epydoc_errors)

        epydoc_output = open('epydoc.out').readlines()
        errors = []
        for line in epydoc_output:
            if line[0] == '|' or line[:2] == '+-':
                errors.append(line)
        if errors:
            fd = open('epydoc.errors', 'w')
            fd.write(''.join(errors))
            fd.close()
            if 1 and email is not None:
                x = os.system(
                    'mail -s "ASE: EpyDoc errors" %s < epydoc.errors' % email)
                assert x == 0

        os.chdir('doc')
        os.mkdir('_build')
        if os.system('PYTHONPATH=%s/lib/python sphinx-build . _build' %
                     tmpdir) != 0:
            raise RuntimeError('Sphinx failed!')
        os.system('cd _build; cp _static/searchtools.js .; ' +
                  'sed -i s/snapshot.tar/%s.tar/ download.html' % version)

        if 1:
            if os.system('PYTHONPATH=%s/lib/python ' % tmpdir +
                         'sphinx-build -b latex . _build 2> error') != 0:
                raise RuntimeError('Sphinx failed!')
            os.system('grep -v "WARNING: unusable reference target found" ' +
                      'error 1>&2')
            # Don't serve mixed content:
            assert os.system(
                'cd _build; find -name "*.html" | xargs sed -i ' +
                '"s%http://cdn.mathjax.org%https://cdn.mathjax.org%"') == 0
            
            os.chdir('_build')
            # os.system('cd ../..; ln -s doc/_static')
            if os.system('make ase-manual.pdf 2>&1') != 0:
                raise RuntimeError('pdflatex failed!')
        else:
            os.chdir('_build')

        assert os.system('mv ../../html epydoc;' +
                         'mv ../../dist/python-ase-%s.tar.gz .' % version) == 0
    
if doc and tarfiledir is not None:
    try:
        os.remove(tarfiledir + '/ase-webpages.tar.gz')
    except OSError:
        pass

build(email)
    
if doc and tarfiledir is not None:
    os.system('cd ..; tar czf %s/ase-webpages.tar.gz _build' % tarfiledir)

os.system('cd; rm -rf ' + tmpdir)
