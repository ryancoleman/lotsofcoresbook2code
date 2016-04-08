.. _hermit:

============================
hermit.hww.de  (Cray XE6) 
============================

Here you find information about the the system
`<http://www.hlrs.de/systems/platforms/cray-xe6-hermit/>`_.

.. note::
   These instructions are up-to-date as of November 18th 2013.

Scalable Python
===============

As the Hermit system is intedend for simulations with thousands of CPU cores,
a special Python interpreter is used here. The scalable Python reduces the 
import time by performing all import related I/O operations with a single CPU 
core and uses MPI for broadcasting the data. As a limitation, all the MPI tasks
have to perform the same **import** statements.

As HLRS does not allow general internet access on compute system, e.g. version 
control repositories cannot be accessed directly (it is possible to setup 
ssh tunnel for some services). Here, we download the scalable Python first to 
to a local machine and use then scp for copying it to Hermit::

  git clone git@gitorious.org:scalable-python/scalable-python.git scalable-python-src
  scp -r scalable-python-src username@hermit1.hww.de:

We will build scalable Python with GNU compilers (other compilers can be used 
for actual GPAW build), so start by changing the default programming 
environment on Hermit::

  module swap PrgEnv-cray PrgEnv-gnu

Due to cross-compile environment in Cray XE6, a normal Python interpreter is 
build for the front-end nodes and the MPI-enabled one for the compute nodes.
The build can be accomplished by the following ``build_gcc`` script

.. literalinclude:: build_gcc.sh

Python packages can now be built on the front-end node with
``/some_path/scalable-python-gcc/bin/python``.

NumPy
=====
As the performance of the HOME filesystem is not very good, we install all the
other components than the pure Python to a disk within the workspace mechanism
of HLRS (with disadvantage that the workspaces expire and have to be 
manually reallocated). Otherwise, no special tricks are needed for installing
NumPy::

  /some_path/scalable-python-gcc/bin/python setup.py install --home=/path_in_workspace


GPAW
====
On Hermit, Intel compiler together with ACML library seemed to give best 
performance for GPAW, in addition HDF5 will be used for parallel I/O. Thus,
load the followgin modules::

  module swap PrgEnv-gnu PrgEnv-intel
  module load acml
  module load hdf5-parallel

The compilation is relatively straightforward, however, as we build NumPy for 
compute nodes it does not work in front-end, and one has to specify NumPy 
include dirs in ``customize.py`` and provide ``--ignore-numpy`` flag when 
building. The system NumPy headers seem to work fine, but safer option is to 
use headers of own NumPy installation

.. literalinclude:: customize_hermit.py

Buid with::

  /some_path/scalable-python-gcc/bin/python setup.py install --home=/path_in_workspace --ignore-numpy

Creating a module
=================
Users can define their own modules for making it easier to setup environment variables. First, create (or edit) a file ``.modulerc`` in the home directory::

  #%Module1.0
  ##

  append-path MODULEPATH $HOME/modules

Now, custom modules can be put to the ``modules`` directory, e.g. GPAW module 
in file ``modules/gpaw``::

  #%Module1.0

  set     prog_root       "/root_folder_of_gpaw_installation"

  # Change PrgEnv to intel
  set confl_prgenvs [list "PrgEnv-cray" "PrgEnv-gnu"]
  foreach prgenv $confl_prgenvs {
    if { [ is-loaded $prgenv] } {
      module swap $prgenv PrgEnv-intel
    }
  }

  setenv GPAW_SETUP_PATH                  "$prog_root/gpaw-setups-0.8.7929"
  prepend-path    PYTHONPATH              "$prog_root/lib/python"
  prepend-path    PATH                    "$prog_root/bin/"

Now, GPAW paths can be set as::

  module load gpaw
