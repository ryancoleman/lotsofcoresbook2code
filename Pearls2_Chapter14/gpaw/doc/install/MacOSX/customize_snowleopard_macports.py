# define_macros += [('GPAW_NO_UNDERSCORE_BLAS', '1')]
# define_macros += [('GPAW_NO_UNDERSCORE_LAPACK', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
# define_macros += [('GPAW_NO_UNDERSCORE_BLACS', '1')]
# define_macros += [('GPAW_NO_UNDERSCORE_SCALAPACK', '1')]
define_macros += [("GPAW_ASYNC",1)]
define_macros += [("GPAW_MPI2",1)]

scalapack = True

libraries = [
           'scalapack',
           'blacsC',
           'blacs',
           'reflapack', # uncomment these two lines 
           'refblas',   # if using vecLib
           'gfortran',
           ]

library_dirs = [
           '/Users/naromero/ScaLAPACK-1.8.0_macportmpich2/lib',
           '/opt/local/lib/gcc45'
           ]

include_dirs += [
    '/opt/local/pkgs/mpich2/include'
    ]

compiler = "/opt/local/bin/mpicc"
mpicompiler = "/opt/local/bin/mpicc"
mpilinker   = "/opt/local/bin/mpicc"
