.. _supernova:

=========
supernova
=========

The supernova machine is a cluster of dual-core AMD Athlon 64 X2 5000+ CPUs,
2.6 GHz processors with 1 GB of memory per core.

Instructions assume **tcsh**, installation under *${HOME}/opt*.
Build the unoptimized numpy/scipy::

  mkdir ${HOME}/opt
  cd ${HOME}/opt

  mkdir -p ${HOME}/opt/python/lib/python2.4/site-packages
  setenv PYTHONPATH ${HOME}/opt/python/lib/python2.4/site-packages

  wget http://dfn.dl.sourceforge.net/sourceforge/numpy/numpy-1.3.0.tar.gz
  wget http://dfn.dl.sourceforge.net/sourceforge/scipy/scipy-0.7.0.tar.gz
  wget http://python-nose.googlecode.com/files/nose-0.11.0.tar.gz
  tar zxf nose-0.11.0.tar.gz
  tar zxf numpy-1.3.0.tar.gz
  tar zxf scipy-0.7.0.tar.gz
  cd nose-0.11.0
  python setup.py install --prefix=${HOME}/opt/python | tee install.log
  cd ../numpy-1.3.0
  python setup.py install --prefix=${HOME}/opt/python | tee install.log
  cd ..
  python -c "import numpy; numpy.test()"

  wget http://www.netlib.org/blas/blas.tgz
  tar zxf blas.tgz
  export BLAS_SRC=${HOME}/opt/BLAS
  wget http://www.netlib.org/lapack/lapack.tgz
  tar zxf lapack.tgz
  export LAPACK_SRC=${HOME}/opt/lapack-3.2.1

  cd scipy-0.7.0
  python setup.py config_fc --fcompiler=gfortran install --prefix=${HOME}/opt/python | tee install.log
  cd ..
  python -c "import scipy; scipy.test()"

Make sure that you have the right mpicc::

  which mpicc
 /usr/local/ompi-1.2.5-pgi/bin/mpicc

and build GPAW (``python setup.py build_ext | tee build_ext.log``)
with this :file:`customize.py` file
(**Note**: instructions valid from the **5232** release)::

  scalapack = True

  compiler = 'gcc'

  extra_compile_args += [
      '-O3',
      '-funroll-all-loops',
      '-fPIC',
      ]

  libraries= []

  mkl_lib_path = '/usr/local/intel/mkl/10.0.011/lib/em64t/'

  library_dirs = [mkl_lib_path]

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


**Note**: is case of problems similar to those found on :ref:`akka` static linking is required.

A gpaw script :file:`test.py` can be submitted like this::

  qsub -l nodes=1:ppn=8 -l walltime=00:30:00 -m abe run.sh

where :file:`run.sh` looks like this::

  #!/bin/sh

  #PBS -m ae
  #PBS -M email@email.com
  #PBS -q long
  #PBS -r n
  #PBS -l nodes=1:ppn=2

  cd $PBS_O_WORKDIR
  echo Running on host `hostname` in directory `pwd`
  NPROCS=`wc -l < $PBS_NODEFILE`
  echo This jobs runs on the following $NPROCS processors:
  cat $PBS_NODEFILE

  export PYTHONPATH=~/opt/gpaw-0.7.5232:~/opt/python-ase-3.1.0.846:${PYTHONPATH}
  export PYTHONPATH=~/opt/python/lib/python2.4/site-packages:${PYTHONPATH}
  export PATH=~/opt/gpaw-0.7.5232/build/bin.linux-x86_64-2.4:${PATH}
  export GPAW_SETUP_PATH=~/opt/gpaw-setups-0.5.3574
  export OMP_NUM_THREADS=1

  mpiexec gpaw-python test.py

Please make sure that your jobs do not run multi-threaded, e.g. for a
job running on ``star237`` do from a login node::

  ssh star237 ps -fL

you should see **1** in the **NLWP** column. Numbers higher then **1**
mean multi-threaded job.

It's convenient to customize as described on the :ref:`parallel_runs` page.
