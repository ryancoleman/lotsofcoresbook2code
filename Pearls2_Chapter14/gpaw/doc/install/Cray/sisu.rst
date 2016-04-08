.. _sisu:

============================
sisu.csc.fi  (Cray XC30) 
============================


.. note::
   These instructions are up-to-date as of July 2014.

GPAW
====

These instructions for GPAW installation use the Scalable Python 
interpreter which reduces drastically the import time in massively parallel
calculations. See the end of this document for installation instructions for 
Scalable Python.

First, the serial version of code is built with serial HDF5 library, e.g. 
for analysis purposes in the front-end::

  module load scalable-python
  module load cray-hdf5

GPAW can be build with a minimal ``customize.py`` (edit the correct paths for
libxc)

.. literalinclude:: customize_sisu.py

Then build the code (no installation at this point yet) with the setup.py 
script::

  python setup.py build_ext

The build of parallel gpaw-python interpreter fails at this point with an error
like::

  build/temp.linux-x86_64-2.6/c/hdf5.o: In function 'h5p_set_fapl_mpio':
  hdf5.c:(.text+0x1bad): undefined reference to 'H5Pset_fapl_mpio'

Next, switch to the parallel version of HDF5 library and do build and install:: 

  module switch cray-hdf5 cray-hdf5-parallel
  python setup.py install --home=path_to_install_prefix

LibXC
=====

Download libxc::

   wget http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.2.0.tar.gz

Configure and make (use GNU environment)::

   ./configure --prefix=install_prefix CC=cc CFLAGS=-fPIC
   make
   make install


Scalable Python
===============

Standard Python interpreter has serious bottleneck in large scale parallel
calculations when many MPI tasks perform the same I/O operations during the
import statetements. Scalable Python `<https://gitorious.org/scalable-python>`_
reduces the import time by having only single MPI task to perform import 
related I/O operations and using then MPI for broadcasting the data.

First, download the source code and switch to GNU environment::

  git clone git@gitorious.org:scalable-python/scalable-python.git
  module switch PrgEnv-cray PrgEnv-gnu

Use the following build script (change the installation prefix to a
proper one):

.. literalinclude:: build_scalable_python.sh

Add then ``install_prefix/bin`` to your PATH, and download and install NumPy::

   export PATH=install_prefix/bin:$PATH
   wget http://sourceforge.net/projects/numpy/files/NumPy/1.8.1/numpy-1.8.1.tar.gz
   tar xzf numpy-1.8.1.tar.gz
   cd numpy-1.8.1
   python setup.py install

