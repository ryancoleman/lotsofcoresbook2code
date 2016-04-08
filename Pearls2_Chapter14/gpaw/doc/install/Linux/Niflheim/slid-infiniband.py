scalapack = True
extra_link_args += ['-cc=gcc']
extra_compile_args += [
  '-cc=gcc',
  '-O2',
  '-m64',
]
libraries = [
  'pathfortran',
  'gfortran',
  'mpiblacsCinit',
  'acml',
  'mpiblacs',
  'scalapack'
  ]
library_dirs = [
  '/opt/pathscale/lib/2.5',
  '/opt/acml-4.0.1/gfortran64/lib',
  '/usr/local/blacs-1.1-24.6.infiniband/lib64',
  '/usr/local/scalapack-1.8.0-1.infiniband/lib64',
  '/usr/local/infinipath-2.0/lib64'
  ]
include_dirs += [
  '/usr/local/infinipath-2.0/include'
 ]
extra_link_args += [
  '-Wl,-rpath=/opt/pathscale/lib/2.5',
  '-Wl,-rpath=/opt/acml-4.0.1/gfortran64/lib',
  '-Wl,-rpath=/usr/local/blacs-1.1-24.6.infiniband/lib64',
  '-Wl,-rpath=/usr/local/scalapack-1.8.0-1.infiniband/lib64',
  '-Wl,-rpath=/usr/local/infinipath-2.0/lib64'
]
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [('GPAW_SECOND_UNDERSCORE_SL_INIT', '1')]

mpicompiler = '/usr/local/infinipath-2.0/bin/mpicc'
mpilinker = mpicompiler
platform_id = 'infiniband'
