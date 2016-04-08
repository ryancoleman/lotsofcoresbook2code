import os

scalapack = True
hdf5 = True
mpicompiler = 'mpicc'
libraries = []

# MKL flags
mkl_flags = os.environ['MKL_SCA_LIBS']
extra_link_args = [mkl_flags]

# HDF5 flags
include_dirs += [os.environ['HDF5_INC_DIR']]
libraries += ['hdf5']
library_dirs += [os.environ['HDF5_LIB_DIR']]

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [("GPAW_ASYNC",1)]
define_macros += [("GPAW_MPI2",1)]

