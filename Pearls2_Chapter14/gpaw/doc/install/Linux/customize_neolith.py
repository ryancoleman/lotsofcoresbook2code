# 19 Nov 2009: problems with scalapack
#scalapack = True

extra_compile_args = ['-fast', '-std=c99', '-fPIC']

compiler = 'icc'

libraries = ['mkl_core', 'mkl_sequential', 'mkl_gf_lp64', 'iomp5']
#libraries = ['mkl_core', 'mkl_intel_thread', 'mkl_gf_lp64', 'mkl_blacs_intelmpi_lp64', 'mkl_scalapack_lp64', 'iomp5']

mkl_lib_path = '/software/intel/mkl/10.2.1.017/lib/em64t/'

library_dirs = [mkl_lib_path]

extra_link_args =['-Wl,-rpath='+mkl_lib_path]

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
