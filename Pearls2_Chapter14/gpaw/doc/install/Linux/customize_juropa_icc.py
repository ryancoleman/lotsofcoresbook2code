compiler = './icc.py'
mpicompiler = './icc.py'
mpilinker = 'MPICH_CC=gcc mpicc'

scalapack = True

library_dirs += ['/opt/intel/Compiler/11.0/074/mkl/lib/em64t']
libraries = ['mkl_intel_lp64' ,'mkl_sequential' ,'mkl_core',
             'mkl_lapack',
             'mkl_scalapack_lp64', 'mkl_blacs_intelmpi_lp64',
             'pthread'
             ]

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [("GPAW_ASYNC",1)]
define_macros += [("GPAW_MPI2",1)]

