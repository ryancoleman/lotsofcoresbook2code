scalapack = False

extra_compile_args = ['-O3', '-std=c99', '-fpic']

compiler = 'gcc'
mpicompiler = '/home/firegam/CAMd/openmpi-1.4.3-1/bin/mpicc'
mpilinker = mpicompiler

libraries = ['mkl_lapack', 'mkl_core', 'mkl_sequential', 'mkl_gf_lp64', 'iomp5']

mkl_lib_path = '/opt/intel/Compiler/11.1/072/mkl/lib/em64t/'
ompi_lib_path = '/home/firegam/CAMd/openmpi-1.4.3-1/lib'

library_dirs = [mkl_lib_path, ompi_lib_path]

extra_link_args =['-Wl,-rpath='+mkl_lib_path+',-rpath='+ompi_lib_path]

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
