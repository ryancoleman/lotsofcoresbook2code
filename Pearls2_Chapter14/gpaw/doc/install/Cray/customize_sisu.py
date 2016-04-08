extra_compile_args = ['-std=c99', '-O3']
compiler = 'cc'
mpicompiler = 'cc'
mpilinker= 'cc'
# edit library and include paths for libxc
include_dirs += ['/homeappl/home/jenkovaa/libxc/sisu/gcc/include']
library_dirs = ['/homeappl/home/jenkovaa/libxc/sisu/gcc/lib']


libraries = ['z', 'xc']

scalapack = True
hdf5 = True

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [("GPAW_ASYNC",1)]
define_macros += [("GPAW_MPI2",1)]
