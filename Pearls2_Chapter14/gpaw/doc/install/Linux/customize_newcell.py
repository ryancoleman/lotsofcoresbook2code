scalapack = False

extra_compile_args = ['-O3', '-std=c99', '-fpic']

compiler = '/afs/crc.nd.edu/x86_64_linux/openmpi/1.3.2/gnu/bin/mpicc'

libraries = ['mkl_def', 'mkl_lapack', 'mkl_core', 'mkl_sequential', 'mkl_gf_lp64', 'iomp5']

mkl_lib_path = '/opt/crc/scilib/mkl/10.1.0.015/lib/em64t/'
ompi_lib_path = '/afs/crc.nd.edu/x86_64_linux/openmpi/1.3.2/gnu/lib/'

library_dirs = [mkl_lib_path, ompi_lib_path]

extra_link_args =['-Wl,-rpath='+mkl_lib_path+',-rpath='+ompi_lib_path]

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
