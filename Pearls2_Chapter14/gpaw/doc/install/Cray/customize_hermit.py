extra_compile_args = ['-std=c99', '-O3'] #, '-O0']
compiler = 'cc'
mpicompiler = 'cc'
mpilinker= 'cc'
libraries = ['acml']
extra_link_args += ['-dynamic']
include_dirs += ['/usr/lib64/python2.6/site-packages/numpy/core/include']
# include_dirs += ['/path_in_workspace/lib/python/numpy/core/include/']

scalapack = True
hdf5 = True

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [("GPAW_ASYNC",1)]
define_macros += [("GPAW_MPI2",1)]
define_macros += [("GPAW_PERFORMANCE_REPORT",1)]
