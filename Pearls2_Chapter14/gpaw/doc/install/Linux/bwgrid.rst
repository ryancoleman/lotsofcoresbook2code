======
bwgrid
======

The `BWgrid <http://www.bw-grid.de/>`__
is an grid of machines located in Baden-Württemberg, Germany.
The installation in Freiburg is a cluster containing 139 dual socket,
quad-core Intel Xenon E5440 CPUs, 2.83GHz processors with 2 GB of memory
per core, 16 dual socket, quad-core Intel Xenon X5550 CPUs, 2.67GHz processors
with 3 GB of memory per core and eight dual socket, six-core Intel Xenon
X5650 CPUs, 2.66GHz processors with 2 GB of memory per core. For more
information visit `<http://www.bfg.uni-freiburg.de/doc/bwgrid>`_.

Building GPAW with Intel compiler
=================================

Use the compiler wrapper file :svn:`~doc/install/Linux/icc.py`

.. literalinclude:: icc.py

Instructions assume **bash**, installation under $HOME/opt.
Load the necessary modules::

  module load devel/python/2.7.2
  module load compiler/intel/12.0
  module load mpi/impi/4.0.2-intel-12.0
  module load numlib/mkl/10.3.5
  module load numlib/python_numpy/1.6.1-python-2.7.2

Internal libxc
--------------

Before revision 10429 libxc was internal. The
:file:`customize.py` had to be changed to
:svn:`~doc/install/Linux/customize_bwgrid_icc.py`

.. literalinclude:: customize_bwgrid_icc.py

External libxc
--------------

After svn revision 10429 libxc has to be included as external library
(see also the `libxc web site <http://www.tddft.org/programs/octopus/wiki/index.php/Libxc:download>`__). To install libxc we assume that MYLIBXCDIR is set to 
the directory where you want to install::

 $ module load compiler/intel/12.0
 $ cd $MYLIBXCDIR
 $ wget http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.0.2.tar.gz
 $ tar -xzvf libxc-2.0.2.tar.gz
 $ cd libxc-2.0.2/
 $ mkdir install
 $ ./configure CFLAGS="-fPIC" --prefix=$PWD/install -enable-shared
 $ make |tee make.log
 $ make install

This will have installed the libs $MYLIBXCDIR/libxc-2.0.2/install/lib 
and the C header
files to $MYLIBXCDIR/libxc-2.0.2/install/include.

We have to modify the file :file:`customize.py` to
:svn:`~doc/install/Linux/customize_bwgrid_icc_libxc.py`

.. literalinclude:: customize_bwgrid_icc_libxc.py

Note that the location of the external libxc on runtime has to be enabled
by setting::

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MYLIBXCDIR/libxc-2.0.2/install/lib

Building and running GPAW
-------------------------

To build GPAW use::

  python setup.py build_ext 2>&1 | tee build_ext.log

and ignore some intermediate warnings.

A gpaw script :file:`test.py` can be submitted to run on 8 cpus like this::

  > gpaw-runscript test.py 8
  using bwg
  run.bwg written
  > qsub run.bwg

