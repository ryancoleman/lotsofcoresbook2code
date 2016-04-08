.. _juropa:

====================================================
juropa.fz-juelich.de   (Intel Xeon, Infiniband, MKL)
====================================================

Here you find information about the the system
http://www.fz-juelich.de/jsc/juropa.

Numpy is installed system wide, so separate installation is not needed.

Building GPAW with gcc
======================

Build GPAW using **gcc** with the configuration file
:svn:`~doc/install/Linux/customize_juropa_gcc.py`.

.. literalinclude:: customize_juropa_gcc.py

and by executing::

  module unload parastation/intel
  module load parastation/gcc

  python setup.py install --prefix='' --home=MY_INSTALLATION_DIR

Building GPAW with Intel compiler
=================================

Use the compiler wrapper file :svn:`~doc/install/Linux/icc.py`

.. literalinclude:: icc.py

Internal libxc
--------------

Before revision 10429 libxc was internal,  
the corresponding 
configuration file is :svn:`~doc/install/Linux/customize_juropa_icc.py`.

.. literalinclude:: customize_juropa_icc.py

External libxc
--------------

After svn revision 10429 libxc has to be included as external library
(see also the `libxc web site <http://www.tddft.org/programs/octopus/wiki/index.php/Libxc:download>`__). To install libxc we assume that MYLIBXCDIR is set to 
the directory where you want to install::

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
:svn:`~doc/install/Linux/customize_juropa_icc_libxc.py`

.. literalinclude:: customize_juropa_icc_libxc.py

Note that the location of the external libxc on runtime has to be enabled
by setting::

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MYLIBXCDIR/libxc-2.0.2/install/lib

Compiling
---------

Now, default parastation/intel module is used so execute only::

  python setup.py install --prefix='' --home=MY_INSTALLATION_DIR

Execution
=========

General execution instructions can be found at
http://www.fz-juelich.de/jsc/juropa/usage/quick-intro.

Example batch job script for GPAW (512 cores, 30 minutes)::

  #!/bin/bash -x
  #MSUB -l nodes=64:ppn=8
  #MSUB -l walltime=0:30:00
  
  cd $PBS_O_WORKDIR
  export PYTHONPATH="MY_INSTALLATION_DIR/ase/lib64/python"
  export PYTHONPATH="$PYTHONPATH":"MY_INSTALLATION_DIR/gpaw/svn/lib64/python"
  export GPAW_SETUP_PATH=SETUP_DIR/gpaw-setups-0.5.3574
  export GPAW_PYTHON=MY_INSTALLATION_DIR/bin/gpaw-python

  export PSP_ONDEMAND=1

  mpiexec -np 512 -x $GPAW_PYTHON my_input.py --sl_default=4,4,64

Note that **-x** flag for *mpiexec* is needed for exporting the environment 
variables to MPI tasks. The environment variable ``PSP_ONDEMAND`` can decrease 
the running time with almost a factor of two with large process counts!

Job scripts can be written also using::

  gpaw-runscript -h

Simultaneous Multi-Threading
============================

SMT_ can be used
to virtually double the number of nodes. A test case did not show
an improvement in performance though.

.. _SMT: http://www2.fz-juelich.de/jsc/juropa/usage/smt

====== ===== === =========
#cores t[s]  SMT date
====== ===== === =========
64     2484  no  9.5.2011
64     2438  no  16.5.2011
128    1081  no  16.5.2011
64     4812  yes 16.5.2011
128    2077  yes 16.5.2011
====== ===== === =========

SMT can be switched on in *gpaw-runscript* via::

  gpaw-runscript -s
