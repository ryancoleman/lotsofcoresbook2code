scalapack = True

libraries = ['scalapack-lam',
             'blacsCinit-lam',
             'blacsF77init-lam',
             'blacs-lam',
             'lapack']

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
