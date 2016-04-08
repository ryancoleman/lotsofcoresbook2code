.. _davinci:

======================
davinci.ssci.liv.ac.uk
======================

The machine is a cluster of dual-core Intel Xeon CPUs, 3.2 GHz
processors with 2 GB of memory per core.

To build (``python setup.py install --home=~/numpy-1.1.0-1``) numpy-1.1.0 add this line to :file:`site.cfg`::

  [DEFAULT]
  library_dirs = /usr/local/Cluster-Apps/intel_mkl_7.0.1.006/mkl701/lib/32

and build GPAW (``PYTHONPATH=${HOME}/dulak/numpy-1.1.0-1/usr/local/lib/python2.5/site-packages python setup.py build_ext``) with this
``customize.py`` file::

  home='/home/haiping'

  extra_compile_args += [
      '-O3'
      ]

  libraries = [
    'mkl',
    'mkl_lapack',
    'guide'
    ]

  library_dirs = [
    '/usr/local/Cluster-Apps/intel_mkl_7.0.1.006/mkl701/lib/32'
    ]

  include_dirs += [
    home+'numpy-1.1.0-1/usr/local/lib/python2.5/site-packages/numpy/core/include'
    ]

A gpaw script :file:`test/CH4.py` can be submitted like this::

  qsub submit.sh

where :file:`submit.sh` looks like this::

  #!/bin/bash
  #
  # Script to submit an mpi job

  # ----------------------------
  # Replace these with the name of the executable 
  # and the parameters it needs

  export home=/home/haiping
  export MYAPP=${home}/gpaw-0.4.2063/build/bin.linux-i686-2.5/gpaw-python
  export MYAPP_FLAGS=${home}/gpaw-0.4.2063/test/CH4.py

  export PYTHONPATH="${home}/numpy-1.1.0-1/usr/local/lib/python2.5/site-packages"
  export PYTHONPATH="${PYTHONPATH}:${home}/gpaw-0.4.2063:${home}/python-ase-3.0.0.358"

  # ---------------------------
  # set the name of the job
  #$ -N CH4

  # request 2 slots
  #$ -pe fatmpi 2


  #################################################################
  #################################################################
  # there shouldn't be a need to change anything below this line

  export MPICH_PROCESS_GROUP=no

  # ---------------------------
  # set up the mpich version to use
  # ---------------------------
  # load the module
  if [ -f  /usr/local/Cluster-Apps/Modules/init/bash ]
  then
          .  /usr/local/Cluster-Apps/Modules/init/bash
          module load default-ethernet
  fi


  #----------------------------
  # set up the parameters for qsub
  # ---------------------------

  #  Mail to user at beginning/end/abort/on suspension
  ####$ -m beas
  #  By default, mail is sent to the submitting user 
  #  Use  $ -M username    to direct mail to another userid 

  # Execute the job from the current working directory
  # Job output will appear in this directory
  #$ -cwd
  #   can use -o dirname to redirect stdout 
  #   can use -e dirname to redirect stderr

  #  Export these environment variables
  #$ -v PATH 
  #$ -v MPI_HOME
  #$ -v LD_LIBRARY_PATH
  #$ -v GPAW_SETUP_PATH
  #$ -v PYTHONPATH
  # Gridengine allocates the max number of free slots and sets the
  # variable $NSLOTS.
  echo "Got $NSLOTS slots."

  # Gridengine sets also $TMPDIR and writes to $TMPDIR/machines the
  # corresponding list of nodes. It also generates some special scripts in
  # $TMPDIR. Therefore, the next two lines are practically canonical:
  #
  #
  export PATH=$TMPDIR:$PATH
  #

  echo "Stack size is "`ulimit -S -s`

  # ---------------------------
  # run the job
  # ---------------------------
  date_start=`date +%s`
  $MPI_HOME/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines  $MYAPP $MYAPP_FLAGS
  date_end=`date +%s`
  seconds=$((date_end-date_start))
  minutes=$((seconds/60))
  seconds=$((seconds-60*minutes))
  hours=$((minutes/60))
  minutes=$((minutes-60*hours))
  echo =========================================================   
  echo SGE job: finished   date = `date`   
  echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
  echo ========================================================= 

It's convenient to customize as in :file:`gpaw-qsub.py` which can be
found at :ref:`parallel_runs`
