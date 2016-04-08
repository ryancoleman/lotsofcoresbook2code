define_macros += [('GPAW_NO_UNDERSCORE_BLAS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_LAPACK', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_BLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_SCALAPACK', '1')]
define_macros += [('GPAW_ASYNC', 1)]
define_macros += [('GPAW_MPI2', 1)]
# define_macros += [('GPAW_MR3',1)] # requires developmental ScaLAPACK
# uncomment two lines below for FLOP rate measurement
# define_macros += [('GPAW_HPM',1)]
define_macros += [('GPAW_PERFORMANCE_REPORT',1)]
# define_macros += [('GPAW_MPI_DEBUG', 1)] # debugging
# define_macros += [('GPAW_OMP',1)] # not really working

scalapack = True
hdf5 = True

# If you are using threading, you probably
# need to change the following library:
# xlomp_ser -> xlsmp
#
# DO NOT INTERCHANGE THE ORDER OF LAPACK
# & ESSL, LAPACK SHOULD BE LINKED FIRST.
#
# Goto appears to be much faster for general
# DGEMM operations, particularly those with:
# alpha != 1.0 and beta != 0.0
#
# Goto is hand-tuned assembly, it will most
# likely always be faster than ESSL-4.x.
# NAR: Goto appears to cause core dumps for
# some problems, use at your own risk.
# Disabling the stackground seems to make
# the problem go away, but this is not 
# recommended.
# --env BG_STACKGUARDENABLE=0
#multi threaded
libraries = [
           'scalapack',
           'lapack',
           'esslsmpbg',
           'xlf90_r',
           'xlopt',
           'xl',
           'xlfmath',
           'xlsmp',
            ]

#single threaded
# libraries = [
#            'scalapack',
#            'lapack',
#            'esslbg',
#            'xlf90_r',
#            'xlopt',
#            'xl',
#            'xlfmath',
#            'xlomp_ser',
#             ]

import os
ibmcmp_base = os.environ['ibmcmp_base']
python_base = '/soft/apps/python/scalable-python-2.6.7-cnk-gcc'

library_dirs = [
           '/soft/libraries/alcf/current/xl/LAPACK/lib',
           '/soft/libraries/alcf/current/xl/SCALAPACK/lib',
           '/soft/libraries/essl/5.1.1-0/lib64',
           '%s/xlf/bg/14.1/bglib64' % ibmcmp_base,
           '%s/xlsmp/bg/3.1/bglib64' % ibmcmp_base,
# plain vanilla Python
#           '/bgsys/tools/Python-2.6/lib64',
# scalable Python 2.6.7
           '%s' % python_base,
           '/soft/libraries/unsupported/hdf5-1.8.8/lib/',
           ]

# plain vanilla Python
# include_dirs += [
#    '/soft/apps/python/python-2.6.6-cnk-gcc/bgsys/tools/Python-2.6/lib64/python2.6/site-packages/numpy/core/include'
#    ]

# scalable Python 2.6.7
include_dirs += [
    '%s/lib/python2.6/site-packages/numpy/core/include' % python_base,
    '/soft/libraries/unsupported/hdf5-1.8.8/include/'
    ]

mpi_libraries = [
#    'mpihpm_smp',
    'hdf5',
    'mpich',
    'opa',
    'mpl',
    'pami',
    'SPI',
    'SPI_cnk',
    'stdc++',
#    'bgpm',
    ]

mpi_library_dirs = [
    '/bgsys/drivers/ppcfloor/comm/xl.legacy/lib',
    '/bgsys/drivers/ppcfloor/comm/sys/lib',
    '/bgsys/drivers/ppcfloor/spi/lib',
    '/soft/perftools/hpctw',
    '/soft/perftools/bgpm/lib',
    ]

extra_link_args = ['-Wl,-export-dynamic']
compiler = "./bgq_xlc.py"
mpicompiler = "./bgq_xlc.py"
mpilinker = "./bgq_xlc_linker.py"
