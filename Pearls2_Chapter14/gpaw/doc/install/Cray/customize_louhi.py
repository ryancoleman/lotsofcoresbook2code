#User provided customizations for the gpaw setup

compiler = 'cc'
mpicompiler = 'cc'
mpilinker= 'cc'

extra_compile_args = ['-std=c99']
libraries = []

scalapack = True
hdf5 = True

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [("GPAW_ASYNC",1)]
define_macros += [("GPAW_MPI2",1)]
