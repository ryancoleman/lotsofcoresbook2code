#!/usr/bin/env python

# Copyright (C) 2007  CAMP
# Please see the accompanying LICENSE file for further information.

from distutils.core import setup, Command
from distutils.command.build_py import build_py as _build_py
from glob import glob
from os.path import join

import os
import sys

import shutil

long_description = """\
ASE is a python package providing an open source Atomic Simulation
Environment in the Python language."""


if sys.version_info < (2, 6, 0, 'final', 0):
    raise SystemExit('Python 2.6 or later is required!')

packages = []
for dirname, dirnames, filenames in os.walk('ase'):
        if '__init__.py' in filenames:
            packages.append(dirname.replace('/', '.'))

package_dir = {'ase': 'ase'}

package_data = {'ase': ['lattice/spacegroup/spacegroup.dat']}


class test(Command):
    description = 'build and run test suite; exit code is number of failures'
    user_options = [('calculators=', 'c',
                     'Comma separated list of calculators to test')]
    
    def __init__(self, dist):
        Command.__init__(self, dist)
        self.sub_commands = ['build']

    def initialize_options(self):
        self.calculators = None

    def finalize_options(self):
        pass

    def run(self):
        self.run_command('build')
        buildcmd = self.get_finalized_command('build')
        sys.path.insert(0, buildcmd.build_lib)

        if self.calculators is not None:
            calculators = self.calculators.split(',')
        elif 'ASE_CALCULATORS' in os.environ:
            calculators = os.environ['ASE_CALCULATORS'].split(',')
        else:
            calculators = []
        from ase.test import test as _test
        testdir = '%s/testase-tempfiles' % buildcmd.build_base
        origcwd = os.getcwd()
        if os.path.isdir(testdir):
            shutil.rmtree(testdir)  # clean before running tests!
        os.mkdir(testdir)
        os.chdir(testdir)
        try:
            results = _test(2, calculators, display=False)
            if results.failures or results.errors:
                print >> sys.stderr, 'Test suite failed'
                raise SystemExit(len(results.failures) + len(results.errors))
        finally:
            os.chdir(origcwd)


class build_py(_build_py):
    """Custom distutils command to build translations."""
    def __init__(self, *args, **kwargs):
        _build_py.__init__(self, *args, **kwargs)
        # Keep list of files to appease bdist_rpm.  We have to keep track of
        # all the installed files for no particular reason.
        self.mofiles = []

    def run(self):
        """Compile translation files (requires gettext)."""
        _build_py.run(self)
        msgfmt = 'msgfmt'
        status = os.system(msgfmt + ' -V')
        if status == 0:
            for pofile in glob('ase/gui/po/*/LC_MESSAGES/ag.po'):
                dirname = join(self.build_lib, os.path.dirname(pofile))
                if not os.path.isdir(dirname):
                    os.makedirs(dirname)
                mofile = join(dirname, 'ag.mo')
                status = os.system('%s -cv %s --output-file=%s 2>&1' %
                                   (msgfmt, pofile, mofile))
                assert status == 0, 'msgfmt failed!'
                self.mofiles.append(mofile)

    def get_outputs(self, *args, **kwargs):
        return _build_py.get_outputs(self, *args, **kwargs) + self.mofiles

# Get the current version number:
execfile('ase/svnversion_io.py')  # write ase/svnversion.py and get svnversion
execfile('ase/version.py')        # get version_base
if svnversion and os.name not in ['ce', 'nt']:
    # MSI accepts only version X.X.X
    version = version_base + '.' + svnversion
else:
    version = version_base

scripts = ['tools/ase-gui', 'tools/ase-db', 'tools/ase-info',
           'tools/ase-build', 'tools/ase-run']
# provide bat executables in the tarball and always for Win
if 'sdist' in sys.argv or os.name in ['ce', 'nt']:
    for s in scripts[:]:
        scripts.append(s + '.bat')

# data_files needs (directory, files-in-this-directory) tuples
data_files = []
for dirname, dirnames, filenames in os.walk('doc'):
    if '.svn' not in dirname:  # skip .svn dirs
        fileslist = []
        for filename in filenames:
            fullname = os.path.join(dirname, filename)
            if '.svn' not in fullname:
                fileslist.append(fullname)
        data_files.append(('share/python-ase/' + dirname, fileslist))

setup(name='python-ase',
      version=version,
      description='Atomic Simulation Environment',
      url='https://wiki.fysik.dtu.dk/ase',
      maintainer='ASE-community',
      maintainer_email='ase-developers@listserv.fysik.dtu.dk',
      license='LGPLv2.1+',
      platforms=['linux'],
      packages=packages,
      package_dir=package_dir,
      package_data=package_data,
      scripts=scripts,
      data_files=data_files,
      long_description=long_description,
      cmdclass={'build_py': build_py,
                'test': test})
