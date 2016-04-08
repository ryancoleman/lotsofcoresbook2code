scalapack = True

extra_compile_args += [
    '-O3'
    ]

libraries = [
           'gfortran',
           'lapack_bgp',
           'scalapack',
           'blacs',
           'lapack_bgp',
           'goto',
           'xlf90_r',
           'xlopt',
           'xl',
           'xlfmath',
           'xlsmp'
           ]

library_dirs = [
           '/home/mdulak/blas-lapack-lib',
           '/home/mdulak/blacs-dev',
           '/home/mdulak/SCALAPACK',
           '/opt/ibmcmp/xlf/bg/11.1/bglib',
           '/opt/ibmcmp/xlsmp/bg/1.7/bglib',
           '/bgsys/drivers/ppcfloor/gnu-linux/lib'
           ]

gpfsdir = '/gpfs/fs2/frontend-13/mdulak'
python_site = 'bgsys/drivers/ppcfloor/gnu-linux'

include_dirs += [gpfsdir+'/Numeric-24.2-1/'+python_site+'/include/python2.5',
                 gpfsdir+'/numpy-1.0.4-1.optimized/'+python_site+'/lib/python2.5/site-packages/numpy/core/include']

extra_compile_args += ['-std=c99']

define_macros += [('GPAW_AIX', '1')]
define_macros += [('GPAW_BGP', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_BLAS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_LAPACK', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_BLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_SCALAPACK', '1')]
