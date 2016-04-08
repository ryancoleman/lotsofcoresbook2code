#!/usr/bin/env python

# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

import os
import sys
import distutils
import distutils.util
from distutils.core import setup, Extension
from glob import glob
from os.path import join

from config import (check_packages, get_system_config, get_parallel_config,
                    get_scalapack_config, get_hdf5_config, check_dependencies,
                    write_configuration, build_interpreter, get_config_vars)


# Get the current version number:
try:
    execfile('gpaw/svnversion_io.py')  # write gpaw/svnversion.py and get
                                       # svnversion
except ValueError:
    svnversion = ''
execfile('gpaw/version.py')        # get version_base
if svnversion:
    version = version_base + '.' + svnversion
else:
    version = version_base

long_description = """\
A grid-based real-space Projector Augmented Wave (PAW) method Density
Functional Theory (DFT) code featuring: Flexible boundary conditions,
k-points and gradient corrected exchange-correlation functionals."""

msg = [' ']

libraries = []
library_dirs = []
include_dirs = []
extra_link_args = []
extra_compile_args = []
runtime_library_dirs = []
extra_objects = []
define_macros = [('NPY_NO_DEPRECATED_API', 7)]
undef_macros = []

mpi_libraries = []
mpi_library_dirs = []
mpi_include_dirs = []
mpi_runtime_library_dirs = []
mpi_define_macros = []

platform_id = ''


packages = []
for dirname, dirnames, filenames in os.walk('gpaw'):
        if '__init__.py' in filenames:
            packages.append(dirname.replace('/', '.'))

include_ase = False
if '--include-ase' in sys.argv:
    include_ase = True
    sys.argv.remove('--include-ase')

import_numpy = True
if '--ignore-numpy' in sys.argv:
    import_numpy = False
    sys.argv.remove('--ignore-numpy')

remove_default_flags = False
if '--remove-default-flags' in sys.argv:
    remove_default_flags = True
    sys.argv.remove('--remove-default-flags')

customize = 'customize.py'
for i, arg in enumerate(sys.argv):
    if arg.startswith('--customize'):
        customize = sys.argv.pop(i).split('=')[1]
        break

check_packages(packages, msg, include_ase, import_numpy)

get_system_config(define_macros, undef_macros,
                  include_dirs, libraries, library_dirs,
                  extra_link_args, extra_compile_args,
                  runtime_library_dirs, extra_objects, msg,
                  import_numpy)

mpicompiler = get_parallel_config(mpi_libraries,
                                  mpi_library_dirs,
                                  mpi_include_dirs,
                                  mpi_runtime_library_dirs,
                                  mpi_define_macros)

mpilinker = mpicompiler
compiler = None

scalapack = False
hdf5 = False

# User provided customizations:
execfile(customize)

if platform_id != '':
    my_platform = distutils.util.get_platform() + '-' + platform_id

    def my_get_platform():
        return my_platform
    distutils.util.get_platform = my_get_platform

if compiler is not None:
    msg += ['* Compiling gpaw with %s' % compiler]
    # A hack to change the used compiler and linker:
    vars = get_config_vars()
    if remove_default_flags:
        for key in ['BASECFLAGS', 'CFLAGS', 'OPT', 'PY_CFLAGS',
                    'CCSHARED', 'CFLAGSFORSHARED', 'LINKFORSHARED',
                    'LIBS', 'SHLIBS']:
            if key in vars:
                value = vars[key].split()
                # remove all gcc flags (causing problems with other compilers)
                for v in list(value):
                    value.remove(v)
                vars[key] = ' '.join(value)
    for key in ['CC', 'LDSHARED']:
        if key in vars:
            value = vars[key].split()
            # first argument is the compiler/linker.  Replace with mpicompiler:
            value[0] = compiler
            vars[key] = ' '.join(value)

custom_interpreter = False
# Check the command line so that custom interpreter is build only with
# 'build', 'build_ext', or 'install':
if mpicompiler is not None:
    for cmd in ['build', 'build_ext', 'install']:
        if cmd in sys.argv:
            custom_interpreter = True
            break

# apply ScaLapack settings
if scalapack:
    get_scalapack_config(define_macros)
    msg.append('* Compiling with ScaLapack')

# distutils clean does not remove the _gpaw.so library and gpaw-python
# binary so do it here:
plat = distutils.util.get_platform()
msg += ['* Architecture: ' + plat]
plat = plat + '-' + sys.version[0:3]
gpawso = 'build/lib.%s/' % plat + '_gpaw.so'
gpawbin = 'build/bin.%s/' % plat + 'gpaw-python'
if 'clean' in sys.argv:
    if os.path.isfile(gpawso):
        print('removing ', gpawso)
        os.remove(gpawso)
    if os.path.isfile(gpawbin):
        print('removing ', gpawbin)
        os.remove(gpawbin)

sources = glob('c/*.c') + ['c/bmgs/bmgs.c']
sources = sources + glob('c/xc/*.c')

check_dependencies(sources)

extension = Extension('_gpaw',
                      sources,
                      libraries=libraries,
                      library_dirs=library_dirs,
                      include_dirs=include_dirs,
                      define_macros=define_macros,
                      undef_macros=undef_macros,
                      extra_link_args=extra_link_args,
                      extra_compile_args=extra_compile_args,
                      runtime_library_dirs=runtime_library_dirs,
                      extra_objects=extra_objects)

extensions = [extension]

if hdf5:
    hdf5_sources = ['c/hdf5.c']
    get_hdf5_config(define_macros)
    msg.append('* Compiling with HDF5')

    hdf5_extension = Extension('_gpaw_hdf5',
                               hdf5_sources,
                               libraries=libraries,
                               library_dirs=library_dirs,
                               include_dirs=include_dirs,
                               define_macros=define_macros,
                               undef_macros=undef_macros,
                               extra_link_args=extra_link_args,
                               extra_compile_args=extra_compile_args,
                               runtime_library_dirs=runtime_library_dirs,
                               extra_objects=extra_objects)
    extensions.append(hdf5_extension)

scripts = [join('tools', script)
           for script in ('gwap', 'gpaw-test', 'gpaw-setup', 'gpaw-basis',
                          'gpaw-mpisim', 'gpaw-runscript',
                          'gpaw-install-setups')]

write_configuration(define_macros, include_dirs, libraries, library_dirs,
                    extra_link_args, extra_compile_args,
                    runtime_library_dirs, extra_objects, mpicompiler,
                    mpi_libraries, mpi_library_dirs, mpi_include_dirs,
                    mpi_runtime_library_dirs, mpi_define_macros)

description = 'An electronic structure code based on the PAW method'

setup(name='gpaw',
      version=version,
      description=description,
      maintainer='GPAW-community',
      maintainer_email='gpaw-developers@listserv.fysik.dtu.dk',
      url='http://wiki.fysik.dtu.dk/gpaw',
      license='GPLv3+',
      platforms=['unix'],
      packages=packages,
      ext_modules=extensions,
      scripts=scripts,
      long_description=long_description)


if custom_interpreter:
    scripts.append('build/bin.%s/' % plat + 'gpaw-python')
    error, par_msg = build_interpreter(
        define_macros, include_dirs, libraries,
        library_dirs, extra_link_args, extra_compile_args,
        runtime_library_dirs, extra_objects,
        mpicompiler, mpilinker, mpi_libraries,
        mpi_library_dirs,
        mpi_include_dirs,
        mpi_runtime_library_dirs, mpi_define_macros)
    msg += par_msg
    # install also gpaw-python
    if 'install' in sys.argv and error == 0:
        setup(name='gpaw',
              version=version,
              description=description,
              maintainer='GPAW-community',
              maintainer_email='gpaw-developers@listserv.fysik.dtu.dk',
              url='http://wiki.fysik.dtu.dk/gpaw',
              license='GPLv3+',
              platforms=['unix'],
              packages=packages,
              ext_modules=[extension],
              scripts=scripts,
              long_description=long_description)

else:
    msg += ['* Only a serial version of gpaw was built!']

# Messages make sense only when building
if 'build' in sys.argv or 'build_ext' in sys.argv or 'install' in sys.argv:
    for line in msg:
        print(line)
