scalapack = True
compiler = 'gcc43'
libraries = [
    'gfortran', 'goto', 'acml',
    'scalapack', 'mpiblacsF77init', 'mpiblacs', 'scalapack',
    # must not link to mpi explicitly: -export-dynamic must be used instead
    ]
library_dirs = [
    '/opt/openmpi/1.3.3-1.el5.fys.gfortran43.4.3.2/lib64',
    '/opt/goto/1.26/1.el5.fys.gfortran43.4.3.2.smp/lib64',
    '/opt/acml/4.3.0/gfortran4364/lib',
    '/opt/blacs/1.1/24.el5.fys.gfortran43.4.3.2.openmpi.1.3.3/lib64',
    '/opt/scalapack/1.8.0/1.el5.fys.gfortran43.4.3.2.openmpi.1.3.3.goto.1.26.acml.4.3.0/lib64'
    ]
include_dirs += ['/opt/openmpi/1.3.3-1.el5.fys.gfortran43.4.3.2/include']
extra_link_args = [
    '-export-dynamic -Wl,-rpath=/opt/openmpi/1.3.3-1.el5.fys.gfortran43.4.3.2/lib64,'
    '-rpath=/opt/goto/1.26/1.el5.fys.gfortran43.4.3.2.smp/lib64,'
    '-rpath=/opt/acml/4.3.0/gfortran4364/lib,'
    '-rpath=/opt/blacs/1.1/24.el5.fys.gfortran43.4.3.2.openmpi.1.3.3/lib64,'
    '-rpath=/opt/scalapack/1.8.0/1.el5.fys.gfortran43.4.3.2.openmpi.1.3.3.goto.1.26.acml.4.3.0/lib64'
    ]
extra_compile_args = ['-O3', '-std=c99', '-funroll-all-loops', '-fPIC']
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]

import tau
import os
tau_path = tau.__file__[0:tau.__file__.find('lib')]
tau_make = tau_path+'lib/Makefile.tau-mpi-pthread-python-pdt'
mpicompiler = "tau_cc.sh -tau_options='-optShared -optCompInst -optVerbose -optMpi' -optTau='-rn Py_RETURN_NONE -i"+os.path.join(os.environ['TAUROOT'], 'include', 'TAU_PYTHON_FIX.h')+"' -tau_makefile="+tau_make
#mpicompiler = "tau_cc.sh -tau_options='-optShared -optCompInst -optVerbose -optMpi' -optTau='-rn Py_RETURN_NONE' -tau_makefile="+tau_make
mpilinker = mpicompiler
compiler = mpicompiler

extra_link_args += ['-Wl,-rpath='+tau_path+'lib/']
platform_id = 'opteron-TAU'
