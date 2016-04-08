.. _akka:

=================
akka.hpc2n.umu.se
=================

The Akka machine (http://www.hpc2n.umu.se/resources/Akka/) is a
cluster of Intel Xeon dual-socket, quad-core L5420 CPUs, 2.5 GHz
processors with 2 GB of memory per core.

On Akka, you need to use the filesystem located under */pfs/* to be
able to write files when running in the queue.  Enable it with
(http://www.hpc2n.umu.se/support/userguide/common/filesystems.html)::

 ln -s /pfs/nobackup$HOME $HOME/pfs

Due to problems with intel mkl
(version *10.0.2.018* gives errors when running on the compute nodes:
*cannot allocate memory for thread-local data: ABORT*)
build numpy using its internal blas/lapack::

 python setup.py install --home=~/pfs/numpy-1.0.4-1

Set these environment variables in the :file:`.bashrc` file::

  export home=~/pfs
  module add openmpi/1.2.6/gcc

  export PYTHONPATH=${home}/gpaw:${home}/ase3k:${home}/numpy-1.0.4-1/lib64/python:
  export GPAW_SETUP_PATH=${home}/gpaw-setups-0.4.2039

  export LD_LIBRARY_PATH=/usr/local/lib
  export PATH=${home}/gpaw/tools:${home}/ase3k/tools:${PATH}

  if [ $PBS_ENVIRONMENT ]; then
        cd $PBS_O_WORKDIR
        export PYTHONPATH=${PBS_O_WORKDIR}:${PYTHONPATH}
        return
  fi

and build GPAW (``python setup.py build_ext``) with this
:file:`customize.py` file (static linking fixes
*cannot allocate memory for thread-local data: ABORT*)::

  scalapack = True

  libraries = []

  extra_compile_args += [
      '-O3'
      ]

  mkl_lib_path = '/usr/local/lib/'

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

A gpaw script :file:`gpaw-script.py` can be submitted like this::

  qsub -l nodes=1:ppn=8 -l walltime=02:00:00 -m abe pbs_submitfile

with the following `pbs_submitfile
<http://www.hpc2n.umu.se/support/userguide/Sarek/src/pbs_submitfile>`_::

  #!/bin/bash
  ###PBS -A SNICXXX-YY-ZZ
  ###PBS -N Parallel
  ###PBS -o test.out
  ###PBS -e test.err
  ###PBS -m ae
  ###PBS -l nodes=2:ppn=8
  ###PBS -l walltime=00:10:00
  ###PBS -l pmem=1900mb # default
  ###PBS -l pvmem=2000mb # default
  
  cd $PBS_O_WORKDIR
  module add openmpi/1.2.6/gcc

  mpiexec ${HOME}/build/bin.linux-x86_64-2.4/gpaw-python gpaw-script.py

It's convenient to customize as described on the :ref:`parallel_runs` page.
