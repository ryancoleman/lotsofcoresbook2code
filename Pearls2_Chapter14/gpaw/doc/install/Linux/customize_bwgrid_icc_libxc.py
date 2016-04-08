compiler = './icc.py'
mpicompiler = './icc.py'
mpilinker = 'MPICH_CC=gcc mpicc'

scalapack = True

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [("GPAW_ASYNC",1)]
define_macros += [("GPAW_MPI2",1)]

libraries += ['mkl_intel_lp64' ,'mkl_sequential' ,'mkl_core',
        'mkl_lapack95_lp64',
        'mkl_scalapack_lp64', 'mkl_blacs_intelmpi_lp64',
        'pthread', 'xc'
        ]

libraries += ['xc']
# change this to your installation directory
LIBXCDIR='/home/mw767/gridpaw/libxc-2.0.2/install/'
library_dirs += [LIBXCDIR + 'lib']
include_dirs += [LIBXCDIR + 'include']
