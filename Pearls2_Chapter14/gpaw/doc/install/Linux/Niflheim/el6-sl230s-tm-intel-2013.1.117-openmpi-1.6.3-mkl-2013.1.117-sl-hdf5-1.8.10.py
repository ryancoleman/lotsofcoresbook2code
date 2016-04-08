nodetype = 'sl230s'
scalapack = True
compiler = 'icc'
libraries =[
    'mkl_def',
    'mkl_scalapack_lp64', 'mkl_intel_lp64', 'mkl_sequential',
    'mkl_core', 'mkl_blacs_openmpi_lp64',
    'hdf5',
    'xc',
    'mpi',
    ]
library_dirs =[
    '/home/opt/common/intel-compilers-2013.1.117/compiler/lib/intel64',
    '/home/opt/el6/' + nodetype + '/openmpi-1.6.3-' + nodetype + '-tm-intel-2013.1.117-1/lib',
    '/home/opt/common/intel-mkl-2013.1.117/mkl/lib/intel64',
    '/home/opt/common/intel-mkl-2013.1.117/compiler/lib/intel64',
    '/home/opt/el6/' + nodetype + '/hdf5-1.8.10-' + nodetype + '-tm-intel-2013.1.117-openmpi-1.6.3-1/lib',
    '/home/opt/el6/' + nodetype + '/libxc-2.0.1-' + nodetype + '-intel-2013.1.117-1/lib',
    ]
include_dirs +=[
    '/home/opt/el6/' + nodetype + '/openmpi-1.6.3-' + nodetype + '-tm-intel-2013.1.117-1/include',
    '/home/opt/el6/' + nodetype + '/hdf5-1.8.10-' + nodetype + '-tm-intel-2013.1.117-openmpi-1.6.3-1/include',
    '/home/opt/el6/' + nodetype + '/libxc-2.0.1-' + nodetype + '-intel-2013.1.117-1/include',
    ]
extra_link_args =[
    '-Wl,-rpath=/home/opt/common/intel-compilers-2013.1.117/compiler/lib/intel64'
    ',-rpath=/home/opt/el6/' + nodetype + '/openmpi-1.6.3-' + nodetype + '-tm-intel-2013.1.117-1/lib'
    ',-rpath=/home/opt/common/intel-mkl-2013.1.117/mkl/lib/intel64'
    ',-rpath=/home/opt/common/intel-mkl-2013.1.117/compiler/lib/intel64'
    ',-rpath=/home/opt/el6/' + nodetype + '/hdf5-1.8.10-' + nodetype + '-tm-intel-2013.1.117-openmpi-1.6.3-1/lib'
    ',-rpath=/home/opt/el6/' + nodetype + '/libxc-2.0.1-' + nodetype + '-intel-2013.1.117-1/lib'
    ]
extra_compile_args =['-xHOST', '-O3', '-ipo', '-std=c99', '-fPIC', '-Wall']
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
mpicompiler = '/home/opt/el6/' + nodetype + '/openmpi-1.6.3-' + nodetype + '-tm-intel-2013.1.117-1/bin/mpicc'
mpilinker = mpicompiler
platform_id = nodetype
hdf5 = True
