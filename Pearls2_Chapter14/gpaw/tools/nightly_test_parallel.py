#!/usr/bin/python
import os
import sys
import time
import glob
import tempfile

def fail(subject, email=None, filename='/dev/null'):
    if email is not None:
        assert os.system('mail -s "%s" %s < %s' %
                         (subject, email, filename)) == 0
    raise SystemExit

if '--dir' in sys.argv:
    i = sys.argv.index('--dir')
    dir = os.path.abspath(sys.argv[i+1])
else:
    dir = None

if '--email' in sys.argv:
    i = sys.argv.index('--email')
    email = sys.argv[i+1]
else:
    email = None

if '--np' in sys.argv:
    i = sys.argv.index('--np')
    np = int(sys.argv[i+1])
else:
    try:
        import multiprocessing
        np = multiprocessing.cpu_count()
    except (ImportError,NotImplementedError):
        np = 2

tmpdir = tempfile.mkdtemp(prefix='gpaw-parallel-', dir=dir)
os.chdir(tmpdir)

# Checkout a fresh version and install:
if os.system('svn checkout ' +
             'https://svn.fysik.dtu.dk/projects/gpaw/trunk gpaw') != 0:
    fail('Checkout of gpaw failed!')
if os.system('svn export ' +
             'https://svn.fysik.dtu.dk/projects/ase/trunk ase') != 0:
    fail('Checkout of ASE failed!')

os.chdir('gpaw')
if os.system('python setup.py ' +
             'install --home=%s 2>&1 | ' % tmpdir +
             'grep -v "c/libxc/src" | tee install.out') != 0:
    fail('Installation failed!', email, 'install.out')

# gpaw installs under libdir
from distutils.sysconfig import get_config_var
libdir = os.path.split(get_config_var('LIBDIR'))[-1]
# import gpaw from where it was installed
sys.path.insert(0, '%s/%s/python' % (tmpdir, libdir))
# and move ase there
os.system('mv ../ase/ase %s/%s/python' % (tmpdir, libdir))

os.system('wget --no-check-certificate --quiet ' +
          'http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-latest.tar.gz')
os.system('tar xvzf gpaw-setups-latest.tar.gz')
setups = tmpdir + '/gpaw/' + glob.glob('gpaw-setups-[0-9]*')[0]

day = time.localtime()[6]
if '--debug' in sys.argv[1:]:
    args = '--debug'
else:
    args = ''

from gpaw.version import version

# Run test-suite:
print('Run')
if os.system('export PYTHONPATH=%s/%s/python:%s/lib/python:${PYTHONPATH}; ' % \
             (tmpdir, libdir, tmpdir) +
             'export GPAW_SETUP_PATH=%s; ' % setups +
             'export OMP_NUM_THREADS=1; ' +
             'mpiexec -np %d ' % np +
             tmpdir + '/bin/gpaw-python ' +
             'tools/gpaw-test --directory=. %s >& test.out' % args) != 0:
    fail('GPAW %s:  Testsuite crashed!' % str(version), email, 'test.out')

try:
    failed = open('failed-tests.txt').readlines()
except IOError:
    pass
else:
    # Send mail:
    n = len(failed)
    if n == 1:
        subject = 'One failed test: ' + failed[0][:-1]
    else:
        subject = '%d failed tests: %s, %s' % (n,
                                               failed[0][:-1], failed[1][:-1])
        if n > 2:
            subject += ', ...'
    subject = 'GPAW %s: ' % str(version) + subject
    fail(subject, email, 'test.out')

print('Done')
os.system('cd; rm -r ' + tmpdir)
