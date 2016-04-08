.. _sepeli:

=============
sepeli.csc.fi
=============

Here you find information about the the system
`<http://raketti.csc.fi/english/research/Computing_services/computing/servers/sepeli_policy>`_.

The installed subversion in sepeli does not support https-protocol, so
one should use a tar file.

Compile like this::

  # use the following modules and define the right python interpreter
  sepeli ~/gpaw/trunk> use mvapich-gnu64
  mvapich-gnu64 is now in use

  MVAPICH environment set
  MPIHOME=/opt/mvapich//gnu64/

  sepeli ~/gpaw/trunk> use ASE
  Atomic Simulation Environment in use
  [ASE is now in use]
  sepeli ~/gpaw/trunk> alias python 'python-pathscale64'
  sepeli ~/gpaw/trunk> unsetenv CC; unsetenv CFLAGS; unsetenv LDFLAGS
  
On runtime you need the following::

  # make shure, that the right acml library is found
  sepeli> setenv LD_LIBRARY_PATH "/opt/acml/gnu64/lib:${LD_LIBRARY_PATH}"

.. Note::

   The compute nodes have different filesystem than the front end
   node. Especially, :envvar:`HOME` and ``$METAWRK`` are
   mounted only on the frontend, so one should place gpaw on 
   ``$WRKDIR``

A sample job script with mvapich (Infiniband) MPI::

   #$ -cwd
   #$ -pe mvapich-gnu64-4 8
   #$ -S /bin/csh
   setenv PYTHONPATH /path_to_ase/:/path_to_gpaw/
   setenv GPAW_SETUP_PATH /path_to_setups/
   setenv PATH "$PATH":/path_to_gpaw-python/
   mpirun -np 8 gpaw-python input.py

In order to use a preinstalled version of gpaw one give the command
``use gpaw`` which sets all the correct environment variables
(:envvar:`PYTHONPATH`, :envvar:`GPAW_SETUP_PATH`, ...)
