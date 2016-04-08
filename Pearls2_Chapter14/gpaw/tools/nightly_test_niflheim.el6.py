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

day = time.localtime()[6]
if '--debug' in sys.argv[1:]:
    args = '--debug'
    np = 2 ** (1 + (day+1) % 3)
else:
    args = ''
    np = 2 ** (1 + day % 3)

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

# allow overwrite of np
if '--np' in sys.argv:
    i = sys.argv.index('--np')
    np = int(sys.argv[i+1])

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
if os.system('source /home/opt/modulefiles/modulefiles_el6.sh&& ' +
             'export PYTHONDONTWRITEBYTECODE=1&& ' +
             'module load NUMPY/1.7.1-1&& ' +
             'module load intel-compilers&& ' +
             'module load openmpi&& ' +
             'python setup.py --remove-default-flags ' +
             '--customize=doc/install/Linux/Niflheim/' +
             'el6-sl230s-tm-gfortran-openmpi-1.6.3-acml-4.4.0-sl-hdf5-1.8.10.py ' +
             'install --home=%s 2>&1 | ' % tmpdir +
             'grep -v "c/libxc/src"') != 0:
    fail('Installation failed!')

os.system('mv ../ase/ase ../lib64/python')

# import gpaw from where it was installed
sys.path.insert(0, '%s/%s/python' % (tmpdir, 'lib64'))

# this import requires numpy!
from gpaw.version import version

os.system('wget --no-check-certificate --quiet ' +
          'http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-latest.tar.gz')
os.system('tar xvzf gpaw-setups-latest.tar.gz')
setups = tmpdir + '/gpaw/' + glob.glob('gpaw-setups-[0-9]*')[0]
    
# Run test-suite:
print('Run')
if os.system('source /home/opt/modulefiles/modulefiles_el6.sh&& ' +
             'export PYTHONDONTWRITEBYTECODE=1&& ' +
             'module load NUMPY/1.7.1-1&& ' +
             'module load SCIPY/0.12.0-1&& ' +
             'module load SCIENTIFICPYTHON&& ' +
             'module load CMR&& ' +
             'module load intel-compilers&& ' +
             'module load openmpi&& ' +
             # libfftw3.so crashes
             'module unload fftw&& ' +  # all fftw must be unloaded
             'export GPAW_FFTWSO=""&& ' +  # and use numpy fftw
             'export PYTHONPATH=%s/lib64/python:$PYTHONPATH; ' % tmpdir +
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
os.system('cd; rm -rf ' + tmpdir)
