scalapack = True
compiler = 'opencc'
libraries =[
    'fortran', 'ffio', 'mv',
    'scalapack', 'mpiblacsF77init', 'mpiblacs', 'scalapack',
    'goto2' ,'acml', 'acml_mv',
    'xc',
    'mpi', 'mpi_f77'
    ]
library_dirs =[
    '/opt/openmpi/1.3.3-1.el5.IB.fys.open64.4.2.3/lib64',
    '/opt/goto2/2.1.13/1.el5.fys.open64.4.2.3.smp/lib64',
    '/opt/acml/4.4.0/open6464/lib',
    '/opt/open64/4.2.3/lib/gcc-lib/x86_64-open64-linux/4.2.3',
    '/opt/blacs/1.1/24.el5.IB.fys.open64.4.2.3.openmpi.1.3.3/lib64',
    '/opt/scalapack/1.8.0/1.el5.IB.fys.open64.4.2.3.openmpi.1.3.3.goto2.2.1.13.acml.4.4.0/lib64',
    '/usr/lib64'
    ]
include_dirs +=['/opt/openmpi/1.3.3-1.el5.IB.fys.open64.4.2.3/include']
extra_link_args =[
    '-Wl,-rpath=/opt/openmpi/1.3.3-1.el5.IB.fys.open64.4.2.3/lib64,'
    '-rpath=/opt/goto2/2.1.13/1.el5.fys.open64.4.2.3.smp/lib64,'
    '-rpath=/opt/acml/4.4.0/open6464/lib,'
    '-rpath=/opt/open64/4.2.3/lib/gcc-lib/x86_64-open64-linux/4.2.3,'
    '-rpath=/opt/blacs/1.1/24.el5.IB.fys.open64.4.2.3.openmpi.1.3.3/lib64,'
    '-rpath=/opt/scalapack/1.8.0/1.el5.IB.fys.open64.4.2.3.openmpi.1.3.3.goto2.2.1.13.acml.4.4.0/lib64,'
    '-rpath=/usr/lib64'  
    ]
extra_compile_args =['-O3', '-std=c99', '-fPIC', '-Wall']
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
mpicompiler = '/opt/openmpi/1.3.3-1.el5.IB.fys.open64.4.2.3/bin/mpicc'
mpilinker = mpicompiler
platform_id = 'opteron-infiniband'
