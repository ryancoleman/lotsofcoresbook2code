scalapack = True

mklpath ='/global/apps/intel/2013.1/mkl'
omppath ='/global/apps/openmpi/1.6.5/intel/13.1'
lxcpath ='/home/pcje/global/apps/libxc-2.2.1-1'

compiler = 'icc'

libraries = ['xc', 'mpi', 'mkl_scalapack_lp64', 'mkl_lapack95_lp64', 'mkl_intel_lp64', 'mkl_sequential', 'mkl_mc', 'mkl_core', 'mkl_def', 'mkl_intel_thread', 'iomp5']
library_dirs += ['%s/lib' % omppath]
library_dirs += ['%s/lib/intel64' % mklpath]
library_dirs += ['%s/lib' % lxcpath]
include_dirs += ['%s/include' % omppath]
include_dirs += ['%s/include' % mklpath]
include_dirs += ['%s/include' % lxcpath]

extra_link_args += ['%s/lib/intel64/libmkl_blacs_openmpi_lp64.a' % mklpath, '%s/lib/intel64/libmkl_blas95_lp64.a' % mklpath]

extra_compile_args += ['-O3', '-std=c99', '-fPIC', '-Wall']

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]

mpicompiler = 'mpicc'
mpilinker = mpicompiler
