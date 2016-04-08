.. _monolith:

========
monolith
========

Here you find information about the the system
`<http://www.nsc.liu.se/systems/retiredsystems/monolith/>`_.

The Monolith machine is a cluster of 2.4Ghz Xeon processors with 2GB of
memory.  The ScaMPI implementation of MPI has a problem, but we can
use MPICH.

Add these two line to the :file:`.modules` file::

  python/2.3.3
  mkl/9.0p18

The Numeric Python module on the system is way too old, so we build
our own version with this :file:`customize.py` file::

  use_system_lapack = 1
  mkl = '/usr/local/intel/ict/l_ict_p_3.0.023/cmkl/9.0'
  lapack_library_dirs = [mkl + '/lib/32']
  lapack_libraries = ['mkl', 'mkl_lapack', 'g2c']
  use_dotblas = 1
  dotblas_include_dirs = [mkl + '/include']
  dotblas_cblas_header = '<mkl_cblas.h>'

Set these environment variables in the :file:`.bashrc` file::

  export PYTHONPATH=$HOME/campos-ase-2.3.4:$HOME/gpaw:$HOME/lib/python/Numeric
  export GPAW_SETUP_PATH=$HOME/setups
  export LD_LIBRARY_PATH=$MKL_ROOT

and build GPAW (``python setup.py build_ext``) with this
:file:`customize.py` file::

  extra_compile_args += ['-w1']
  mpicompiler = 'icc -Nmpich'
  custom_interpreter = True
  compiler = 'icc'

Jobs can be submitted like this::

  qsub -l nodes=2:ppn=2 -A <account> -l walltime=2:00:00 \
       -m abe run.sh

where :file:`run.sh` looks like this::

  cd $PBS_O_WORKDIR
  mpirun $HOME/gpaw/build/bin.linux-i686-2.3/gpaw-python gpaw-script.py
