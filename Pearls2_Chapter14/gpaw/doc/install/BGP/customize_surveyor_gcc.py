define_macros += [('GPAW_NO_UNDERSCORE_BLAS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_LAPACK', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_BLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_SCALAPACK', '1')]
define_macros += [('GPAW_ASYNC',1)]
define_macros += [('GPAW_MPI2',1)]
# define_macros += [('GPAW_MR3',1)] # requires experimental ScaLAPACK
# uncomment two lines below for FLOP rate measurement
# define_macros += [('GPAW_HPM',1)]
# define_macros += [('GPAW_PERFORMANCE_REPORT',1)]
define_macros += [('GPAW_MPI_DEBUG',1)] # debugging
# define_macros += [('GPAW_OMP',1)] # not really working

hdf5 = True
scalapack = True

# If you are using threading, you probably
# need to change the following library:
# xlomp_ser -> xlsmp
#
# DO NOT INTERCHANGE THE ORDER OF LAPACK
# & ESSL, LAPACK SHOULD BE LINKED FIRST.
#
# Goto BLAS is broken on BG/P. It should not be used.
#

libraries = [
#           'scalapackmr3',
           'scalapack',
           'blacsCinit_MPI-BGP-0',
           'blacs_MPI-BGP-0',
           'lapack_bgp',
           'esslbg',
           'hdf5',
           'xlf90_r',
           'xlopt',
           'xl',
           'xlfmath',
           'xlomp_ser',
#           'hpm',
           ]

#          make sure XL library_dirs below match XL compiler version
#          (e.g. aug2010, jan2011) used in mpilinker variable

library_dirs = [
           '/soft/apps/SCALAPACK-dev-r98',
           '/soft/apps/LAPACK',
           '/soft/apps/ESSL-4.4.1-1/lib',
           '/soft/apps/ibmcmp-aug2011/xlf/bg/11.1/bglib',
           '/soft/apps/ibmcmp-aug2011/xlsmp/bg/1.7/bglib',
           '/bgsys/drivers/ppcfloor/gnu-linux/lib',
#           '/soft/apps/UPC/lib',
           '/soft/apps/hdf5-1.8.0/lib',
           ]

include_dirs += [
    '/soft/apps/python/python-2.6-cnk-gcc/numpy-1.3.0/lib/python2.6/site-packages/numpy/core/include',
    '/soft/apps/hdf5-1.8.0/include'
    ]

# TAU library below needed for automatic instrumentation only

mpi_libraries = [
#   'TAU',
    'mpich.cnk',
    'opa',
    'dcmf.cnk',
    'dcmfcoll.cnk',
    'SPI.cna',
    ]

mpi_library_dirs = [
    '/soft/apps/tau/tau-2.19.2/bgp/lib/bindings-bgptimers-gnu-mpi-python-pdt',
    '/bgsys/drivers/ppcfloor/comm/default/lib',
    '/bgsys/drivers/ppcfloor/comm/sys/lib',
    '/bgsys/drivers/ppcfloor/runtime/SPI',
    ]

extra_link_args += ['-Wl,-export-dynamic'] # make symbols in *.a visible to *.so, needed for TAU
compiler = "bgp_gcc.py"
mpicompiler = "bgp_gcc.py"
mpilinker   = "bgp_gcc_linker.py"
