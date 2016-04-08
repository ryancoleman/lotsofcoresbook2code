#!/usr/bin/python

import os
import sys
import time
import glob
import trace
import tempfile

def send_email(subject, filename='/dev/null'):
    assert os.system('mail -s "%s" %s < %s' %
                     (subject, email, filename)) == 0

def fail(msg, email=None, filename='/dev/null'):
    if email is not None:
        send_email(msg, filename)
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

tmpdir = tempfile.mkdtemp(prefix='gpaw-', dir=dir)
os.chdir(tmpdir)

# Checkout a fresh version and install:
if os.system('svn checkout ' +
             'https://svn.fysik.dtu.dk/projects/gpaw/trunk gpaw') != 0:
    fail('Checkout of gpaw failed!')

if os.system('svn export ' +
             'https://svn.fysik.dtu.dk/projects/ase/trunk ase') != 0:
    fail('Checkout of ASE failed!')

os.chdir('gpaw')

if os.system('python setup.py install --home=%s ' % tmpdir +
             '2>&1 | grep -v "c/libxc/src" | tee install.out') != 0:
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
if day % 2:
    sys.argv.append('--debug')

from gpaw import setup_paths
setup_paths.insert(0, setups)

from gpaw.version import version

# Run test-suite:
from gpaw.test import TestRunner, tests
os.mkdir('gpaw-test')
os.chdir('gpaw-test')
out = open('test.out', 'w')
#tests = ['ase3k.py']
failed = TestRunner(tests, stream=out).run()
out.close()
if failed:
    # Send mail:
    n = len(failed)
    if n == 1:
        subject = 'One failed test: ' + failed[0]
    else:
        subject = ('%d failed tests: %s, %s' %
                   (n, failed[0], failed[1]))
        if n > 2:
            subject += ', ...'
    subject = 'GPAW %s: ' % str(version) + subject
    fail(subject, email, 'test.out')

def count(dir, pattern):
    p = os.popen('wc -l `find %s -name %s` | tail -1' % (dir, pattern), 'r')
    return int(p.read().split()[0])

os.chdir('..')
if 0:  # revision 10429 - libxc merged
    libxc = count('c/libxc', '\\*.[ch]')
else:
    libxc = 0
ch = count('c', '\\*.[ch]') - libxc
test = count('gpaw/test', '\\*.py')
py = count('gpaw', '\\*.py') - test

"""
import pylab
# Update the stat.dat file:
dir = '/tmp/nightly-test/'
f = open(dir + 'stat.dat', 'a')
print >> f, pylab.epoch2num(time.time()), libxc, ch, py, test
f.close()

# Construct the stat.png file:
lines = open(dir + 'stat.dat').readlines()
date, libxc, c, code, test = zip(*[[float(x) for x in line.split()]
                                   for line in lines[1:]])
date = pylab.array(date)
code = pylab.array(code)
test = pylab.array(test)
c = pylab.array(c)

def polygon(x, y1, y2, *args, **kwargs):
    x = pylab.concatenate((x, x[::-1]))
    y = pylab.concatenate((y1, y2[::-1]))
    pylab.fill(x, y, *args, **kwargs)

fig = pylab.figure()
ax = fig.add_subplot(111)
polygon(date, code + test, code + test + c,
        facecolor='r', label='C-code')
polygon(date, code, code + test,
        facecolor='y', label='Tests')
polygon(date, [0] * len(date), code,
        facecolor='g', label='Python-code')
polygon(date, [0] * len(date), [0] * len(date),
        facecolor='b', label='Fortran-code')

months = pylab.MonthLocator()
months3 = pylab.MonthLocator(interval=3)
month_year_fmt = pylab.DateFormatter("%b '%y")

ax.xaxis.set_major_locator(months3)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_major_formatter(month_year_fmt)
labels = ax.get_xticklabels()
pylab.setp(labels, rotation=30)
pylab.axis('tight')
pylab.legend(loc='upper left')
pylab.title('Number of lines')
pylab.savefig(dir + 'stat.png')
"""

print('Done')
os.system('cd; rm -rf ' + tmpdir)
