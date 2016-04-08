.. _murska:

=========================================================
murska.csc.fi  (AMD Opteron, Infiniband, ACML)
=========================================================

Here you find information about the the system
`<http://www.csc.fi/english/research/Computing_services/computing/servers/murska>`_.

Installation of user's packages on murska is recommended under /v/users/$USER/appl/.

We want to use python2.4 and gcc compiler::

  > module load python
  > module swap PrgEnv-pgi PrgEnv-gnu

and use this :file:`customize.py`::

  scalapack = True

  extra_compile_args =['-O3', '-std=c99']

  libraries =['gfortran','acml','scalapack','mpiblacsF77init','mpiblacs','scalapack']
  library_dirs =[
        '/v/linux26_x86_64/opt/blacs/1.1gnu/hpmpi/lib64',
        '/v/linux26_x86_64/opt/scalapack/1.8.0gnu/scalapack-1.8.0'
        ]

  define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
  define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]

Then, compile GPAW ``python setup.py build_ext``.

A sample job script::

  #!/bin/csh

  #BSUB -n 4
  #BSUB -W 0:10
  #BSUB -J jobname_%J
  #BSUB -e jobname_err_%J
  #BSUB -o jobname_out_%J

  # If you install you personal version of gpaw under /v/users/$USER/appl/
  # load the required modules and set the environment variables PYTHONPATH, etc.
  module load ASE/svn
  module swap PrgEnv-pgi PrgEnv-gnu
  module load gpaw-setups
  setenv PYTHONPATH /v/users/$USER/appl/gpaw:$PYTHONPATH
  setenv PATH /v/users/$USER/appl/gpaw/build/bin.linux-x86_64-2.4:$PATH

  # Alternatively, use a preinstalled version of gpaw load just gpaw/svn module
  # which sets all the correct environment variables
  # PYTHONPATH, PATH, GPAW_SETUP_PATH, etc.
  # module load gpaw/svn

  mpirun -srun gpaw-python input.py

Murska uses LSF-HPC batch system where jobs are submitted as (note the
stdin redirection)::

  > bsub < input.py
