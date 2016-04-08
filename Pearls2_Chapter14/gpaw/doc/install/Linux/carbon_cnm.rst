.. _carbon_cnm:

======
carbon
======

Here you find information about the the system
`<http://www.top500.org/system/9025>`_.

The Carbon machine is a cluster of dual socket, quad-core Intel Xeon
5355 CPUs, 2.66 GHz processors with 2 GB of memory per core.

**Warning**: numpy build instructions have not been tested recently and may not work.
Please use the (unoptimized) system default numpy.

To build (``python setup.py install --home=~/numpy-1.0.4-1``)
numpy-1.0.4 add these lines to :file:`site.cfg`::

  [DEFAULT]
  library_dirs = /usr/local/lib:/opt/intel/mkl/10.0.2.018/lib/em64t
  include_dirs = /usr/local/include:/opt/intel/mkl/10.0.2.018/include

and, in :file:`numpy/distutils/system_info.py` change the line::

  _lib_mkl = ['mkl','vml','guide']

into::

  _lib_mkl = ['mkl','guide']

and the line::

  lapack_libs = self.get_libs('lapack_libs',['mkl_lapack32','mkl_lapack64'])

into::

  lapack_libs = self.get_libs('lapack_libs',['mkl_lapack'])

Set these environment variables in the :file:`.bashrc` file::

  export OMPI_CC=gcc
  export OMP_NUM_THREADS=1

  export PYTHONPATH=${HOME}/gpaw:${HOME}/ase3k:${HOME}/numpy-1.0.4-1/lib64/python:
  export GPAW_SETUP_PATH=${HOME}/gpaw-setups-0.4.2039

  export LD_LIBRARY_PATH=/usr/lib64/openmpi:/opt/intel/mkl/10.0.2.018/lib/em64t
  export PATH=${HOME}/gpaw/tools:${HOME}/ase3k/tools:/usr/share/openmpi/bin64:${PATH}

  if [ $PBS_ENVIRONMENT ]; then
        cd $PBS_O_WORKDIR
        export PYTHONPATH=${PBS_O_WORKDIR}:${PYTHONPATH}
        return
  fi

and build GPAW (``python setup.py build_ext``) with this
:file:`customize.py` file::

  scalapack = True

  extra_compile_args += [
      '-O3'
      ]

  libraries= []

  mkl_lib_path = '/opt/intel/mkl/10.0.2.018/lib/em64t/'

  extra_link_args = [
  mkl_lib_path+'libmkl_intel_lp64.a',
  mkl_lib_path+'libmkl_sequential.a',
  mkl_lib_path+'libmkl_core.a',
  mkl_lib_path+'libmkl_blacs_openmpi_lp64.a',
  mkl_lib_path+'libmkl_scalapack.a',
  mkl_lib_path+'libmkl_blacs_openmpi_lp64.a',
  mkl_lib_path+'libmkl_intel_lp64.a',
  mkl_lib_path+'libmkl_sequential.a',
  mkl_lib_path+'libmkl_core.a',
  mkl_lib_path+'libmkl_intel_lp64.a',
  mkl_lib_path+'libmkl_sequential.a',
  mkl_lib_path+'libmkl_core.a',
  ]

  define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
  define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]

**Note**: due to linking problems similar to those found on :ref:`akka` static linking is required.

A gpaw script :file:`gpaw-script.py` can be submitted like this::

  qsub -l nodes=1:ppn=8 -l walltime=02:00:00 \
       -m abe run.sh

where :file:`run.sh` looks like this::

  cd $PBS_O_WORKDIR
  mpirun -machinefile $PBS_NODEFILE -np 8 -x OMP_NUM_THREADS \
         $HOME/gpaw/build/bin.linux-x86_64-2.4/gpaw-python gpaw-script.py

Please make sure that your jobs do not run multi-threaded, e.g. for a
job running on ``n090`` do from a login node::

  ssh n090 ps -fL

you should see **1** in the **NLWP** column. Numbers higher then **1**
mean multi-threaded job.

It's convenient to customize as described on the :ref:`parallel_runs` page.
