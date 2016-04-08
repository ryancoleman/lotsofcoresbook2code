scalapack = True
extra_compile_args += ['-O3']
libraries = [
  'gfortran',
  'mpiblacsCinit',
  'mpiblacs',
  'acml',
  'scalapack',
  'mpi_f77'
  ]
library_dirs = [
  '/opt/acml-4.0.1/gfortran64/lib',
  '/usr/local/blacs-1.1-24.56.gfortran/lib64',
  '/usr/local/scalapack-1.8.0-1.gfortran.acml/lib64',
  '/usr/local/openmpi-1.2.5-gfortran/lib64'
  ]
include_dirs += [
  '/usr/local/openmpi-1.2.5-gfortran/include'
  ]
extra_link_args += [
  '-Wl,-rpath=/opt/acml-4.0.1/gfortran64/lib',
  '-Wl,-rpath=/usr/local/blacs-1.1-24.56.gfortran/lib64',
  '-Wl,-rpath=/usr/local/scalapack-1.8.0-1.gfortran.acml/lib64',
  '-Wl,-rpath=/usr/local/openmpi-1.2.5-gfortran/lib64'
  ]
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
mpicompiler = '/usr/local/openmpi-1.2.5-gfortran/bin/mpicc'
mpilinker = mpicompiler
platform_id = 'ethernet'
