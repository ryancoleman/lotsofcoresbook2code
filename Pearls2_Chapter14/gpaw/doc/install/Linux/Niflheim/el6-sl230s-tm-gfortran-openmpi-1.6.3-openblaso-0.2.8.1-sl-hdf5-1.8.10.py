nodetype = 'sl230s'
scalapack = True
compiler = 'gcc'
libraries =[
    'gfortran',
    'scalapack',
    'mpiblacs',
    'mpiblacsCinit',
    'openblaso',
    'hdf5',
    'xc',
    'mpi',
    'mpi_f77',
    ]
library_dirs =[
    '/home/opt/el6/' + nodetype + '/openmpi-1.6.3-' + nodetype + '-tm-gfortran-1/lib',
    '/home/opt/el6/' + nodetype + '/blacs-1.1-' + nodetype + '-tm-gfortran-openmpi-1.6.3-1/lib',
    '/home/opt/el6/' + nodetype + '/scalapack-2.0.2-' + nodetype + '-tm-gfortran-openmpi-1.6.3-acml-4.4.0-1/lib',
    '/home/opt/el6/common/openblas-0.2.8-1/lib64',
    '/home/opt/el6/' + nodetype + '/hdf5-1.8.10-' + nodetype + '-tm-gfortran-openmpi-1.6.3-1/lib',
    '/home/opt/el6/' + nodetype + '/libxc-2.2.1-' + nodetype + '-gfortran-1/lib',
    ]
include_dirs +=[
    '/home/opt/el6/' + nodetype + '/openmpi-1.6.3-' + nodetype + '-tm-gfortran-1/include',
    '/home/opt/el6/' + nodetype + '/hdf5-1.8.10-' + nodetype + '-tm-gfortran-openmpi-1.6.3-1/include',
    '/home/opt/el6/' + nodetype + '/libxc-2.2.1-' + nodetype + '-gfortran-1/include',
    ]
extra_link_args =[
    '-Wl,-rpath=/home/opt/el6/' + nodetype + '/openmpi-1.6.3-' + nodetype + '-tm-gfortran-1/lib'
    ',-rpath=/home/opt/el6/' + nodetype + '/blacs-1.1-' + nodetype + '-tm-gfortran-openmpi-1.6.3-1/lib'
    ',-rpath=/home/opt/el6/' + nodetype + '/scalapack-2.0.2-' + nodetype + '-tm-gfortran-openmpi-1.6.3-acml-4.4.0-1/lib'
    ',-rpath=/home/opt/el6/common/openblas-0.2.8-1/lib64'
    ',-rpath=/home/opt/el6/' + nodetype + '/hdf5-1.8.10-' + nodetype + '-tm-gfortran-openmpi-1.6.3-1/lib'
    ',-rpath=/home/opt/el6/' + nodetype + '/libxc-2.2.1-' + nodetype + '-gfortran-1/lib'
    ]
extra_compile_args =['-O3', '-std=c99', '-fPIC', '-Wall']
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
mpicompiler = '/home/opt/el6/' + nodetype + '/openmpi-1.6.3-' + nodetype + '-tm-gfortran-1/bin/mpicc'
mpilinker = mpicompiler
platform_id = nodetype
hdf5 = True
