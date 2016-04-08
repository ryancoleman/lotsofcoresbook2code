scalapack = False

extra_compile_args = ['-O3', '-std=c99', '-fpic']

compiler = 'gcc'
mpicompiler = '/usr/local/openmpi-1.3.3/bin/mpicc'
mpilinker = mpicompiler

libraries = ['mkl_lapack', 'mkl_core', 'mkl_sequential', 'mkl_gf', 'iomp5']

mkl_lib_path = '/opt/intel/mkl/10.2.1.017/lib/32'
ompi_lib_path = '/usr/local/openmpi-1.3.3/lib'

library_dirs = [mkl_lib_path, ompi_lib_path]

extra_link_args =['-Wl,-rpath='+mkl_lib_path+',-rpath='+ompi_lib_path]

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
