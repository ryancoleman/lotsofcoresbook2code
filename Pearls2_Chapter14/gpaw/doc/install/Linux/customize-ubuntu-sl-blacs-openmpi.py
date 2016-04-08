scalapack = True

libraries = ['scalapack-openmpi',
             'blacsCinit-openmpi',
             'blacsF77init-openmpi',
             'blacs-openmpi',
             'xc',
             'blas',
             'lapack']

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
