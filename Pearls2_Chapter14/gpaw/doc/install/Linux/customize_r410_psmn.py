scalapack = False

compiler = 'gcc'

extra_compile_args += [
    '-O3',
    '-funroll-all-loops',
    '-fPIC',
    ]

libraries = ['gfortran']

mpi_prefix = '/softs/openmpi-gnu/'
blas_lib_path = '/home/tjiang/softs/acml-4.0.1/gfortran64/lib/'
lapack_lib_path = blas_lib_path

library_dirs = [mpi_prefix + 'lib']
include_dirs += [mpi_prefix + 'include'] # includes alreay numpy's include

extra_link_args = [
    blas_lib_path+'libacml.a',
    lapack_lib_path+'libacml.a',
    '-Wl,-rpath=' + mpi_prefix + 'lib,'
    '-rpath=' + blas_lib_path
    ]

mpicompiler= mpi_prefix + 'bin/mpicc'
mpilinker = mpicompiler
