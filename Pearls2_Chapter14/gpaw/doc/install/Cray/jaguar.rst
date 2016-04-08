.. _jaguar:

==================
jaguar  (Cray XT5)
==================

Here you find information about the the system
http://www.nccs.gov/computing-resources/jaguar/.

The current operating system in Cray XT4/XT5 compute nodes, Compute Linux
Environment (CLE) has some limitations, most notably it does not
support shared libraries. In order to use python in CLE some
modifications to the standard python are needed. Instructions below
assume **tcsh**.

The installations process of python and numpy can be performed with the
script :svn:`~doc/install/Cray/make_python_numpy`::

  ./make_python_numpy |& tee all.log

whose details are given below. **Note**: One may want to change the
installation paths in the beginning of *make_python_numpy*.

Set the correct C compiler and flags, e.g.::

  module swap PrgEnv-pgi PrgEnv-gnu
  setenv CC cc
  setenv CXX CC
  setenv OPT '-O3 -funroll-all-loops'

The following modules are loaded::

  module avail
   Currently Loaded Modulefiles:
   1) modules/3.1.6
   2) DefApps
   3) torque/2.2.0-snap.200707311754
   4) moab/5.2.4
   5) xtpe-quadcore
   6) MySQL/5.0.45
   7) xt-service/2.1.50HD
   8) xt-libc/2.1.50HD
   9) xt-os/2.1.50HD
  10) xt-boot/2.1.50HD
  11) xt-lustre-ss/2.1.50HD.PS04.lus.1.6.5.steve.8103_1.6.5
  12) xtpe-target-cnl
  13) Base-opts/2.1.50HD
  14) PrgEnv-gnu/2.1.50HD
  15) xt-asyncpe/2.3
  16) xt-pe/2.1.50HD
  17) xt-mpt/3.1.0
  18) xt-libsci/10.3.1
  19) fftw/3.1.1
  20) gcc/4.2.0.quadcore

The recommended place for user's applications is under :envvar:`HOME`::

  cd
  mkdir -p sw/xt5
  cd sw/xt5
  set sw_home=~/sw/xt5
  wget http://www.python.org/ftp/python/2.5.4/Python-2.5.4.tar.bz2
  wget http://sunet.dl.sourceforge.net/sourceforge/expat/expat-2.0.1.tar.gz
  wget http://www.zlib.net/zlib-1.2.3.tar.bz2
  tar jxf Python-2.5.4.tar.bz2
  tar zxf expat-2.0.1.tar.gz
  tar jxf zlib-1.2.3.tar.bz2
  wget http://python-nose.googlecode.com/files/nose-0.11.0.tar.gz
  tar zxf nose-0.11.0.tar.gz
  wget http://dfn.dl.sourceforge.net/sourceforge/numpy/numpy-1.2.1.tar.gz
  tar zxf numpy-1.2.1.tar.gz

Before installing a special python, expat_ and zlib_
which are needed by GPAW,
but which are not included in the python distribution.
The installation is based on instructions from
http://yt.enzotools.org/wiki/CrayXT5Installation.

.. _expat: http://expat.sourceforge.net/
.. _zlib: http://www.zlib.net/  

Install expat::

  cd ${sw_home}
  setenv EXPAT_DIR ${sw_home}/expat-2.0.1-1
  cd expat-2.0.1
  ./configure --disable-shared --prefix=${EXPAT_DIR}
  make
  make install

Install zlib::

  cd ${sw_home}
  setenv ZLIB_DIR ${sw_home}/zlib-1.2.3-1
  cd zlib-1.2.3
  ./configure --prefix=${ZLIB_DIR}
  make # ignore error: /usr/lib/../lib64/libc.a: could not read symbols: Bad value
  make install

Next, one can proceed with the actual python installation. The
following instructions are tested with python 2.5.4:

- enter the python source directory::

   cd ${sw_home}
   setenv PYTHON_DIR ${sw_home}/Python-2.5.4-1
   cd Python-2.5.4

- create a special dynamic loader for correct resolution of namespaces::

   wget --no-check-certificate http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Cray/dynload_redstorm.c -O Python/dynload_jaguar.c

- run :file:`configure`::

   ./configure --prefix=${PYTHON_DIR} SO=.a DYNLOADFILE=dynload_jaguar.o MACHDEP=jaguar --host=x86_64-unknown-linux-gnu --disable-sockets --disable-ssl --enable-static --disable-shared | tee config.log

- in order to use ``distutils`` append the :file:`Lib/distutils/unixccompiler.py` file, so that static libraries are created instead of shared ones::

   wget --no-check-certificate http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Cray/linkforshared.py
   cat Lib/distutils/unixccompiler.py linkforshared.py > unixccompiler.py
   mv unixccompiler.py  Lib/distutils

- specify which modules will be statically linked in to the python interpreter
  by editing :file:`Modules/Setup`::

   mv Modules/Setup Modules/Setup.orig
   wget --no-check-certificate http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Cray/Setup_jaguar -O Modules/Setup
   touch Modules/Setup

  **Note**: sha modules are required by numpy, so the following lines should be present in Modules/Setup::

   _sha shamodule.c
   _sha256 sha256module.c
   _sha512 sha512module.c

-  modify :file:`Lib/locale.py` as described at `<http://yt.enzotools.org/wiki/CrayXT5Installation>`_ (is it really needed?),

- build and install::

   make | tee make.log
   # ignore errors like:
   # *** WARNING: renaming "_ctypes" since importing it failed: dynamic module does not define init function (init_ctypes)
   make install | tee make_install.log   

- build numpy::

   cd ${sw_home}
   cd numpy-1.2.1
   ${PYTHON_DIR}/bin/python setup.py install | tee install.log

  **Note**: numpy 1.3.0 gives::

   # ImportError: No module named select

- append numpy to pythons's :file:`Modules/Setup`::

   cd ${sw_home}/Python-2.5.4
   cat ../numpy-1.2.1/install.log | grep Append | cut -d ":" -f 2 | sed -n 's/ *//p' > append
   cat Modules/Setup append > Setup
   mv Setup Modules

  example output::

   cat append
   multiarray /autofs/na1_home/farberow/sw/xt5/numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/core/multiarray.a
   umath /autofs/na1_home/farberow/sw/xt5/numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/core/umath.a
   _sort /autofs/na1_home/farberow/sw/xt5/numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/core/_sort.a
   scalarmath /autofs/na1_home/farberow/sw/xt5/numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/core/scalarmath.a
   _compiled_base /autofs/na1_home/farberow/sw/xt5/numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/lib/_compiled_base.a
   _capi /autofs/na1_home/farberow/sw/xt5/numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/numarray/_capi.a
   fftpack_lite /autofs/na1_home/farberow/sw/xt5/numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/fft/fftpack_lite.a
   lapack_lite /autofs/na1_home/farberow/sw/xt5/numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/linalg/lapack_lite.a
   mtrand /autofs/na1_home/farberow/sw/xt5/numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/random/mtrand.a

- rebuild python::

   make | tee make2.log
   make install | tee make_install2.log

On jaguar only */tmp/work/$USER* filesystem is available for batch jobs.
**Note**: that this space is cleaning periodically
http://www.nccs.gov/computing-resources/jaguar/file-systems/.
Test python/numpy::

 cp -r ${PYTHON_DIR} /tmp/work/$USER
 cp -r ${sw_home}/nose-0.11.0 /tmp/work/$USER
 cd /tmp/work/$USER

 cat <<EOF > ./numpyTest.py
 import numpy
 from numpy.core.multiarray import dot
 b = numpy.ones(13, numpy.complex)
 d = dot(b, b)
 print 'Hello'
 numpy.test()
 EOF

 cat <<EOF > ./numpyTest.pbs
 #!/bin/bash
 #PBS -l walltime=00:10:00,size=8
 #PBS -N numpyTest
 #PBS -A XXXXXX
 #PBS -j oe

 export PYTHONHOME=/tmp/work/$USER/Python-2.5.4-1
 export PYTHONPATH=/tmp/work/$USER/nose-0.11.0

 cd /tmp/work/$USER
 env | grep PYTHON
 env | grep LD_LIBRARY_PATH
 aprun -n1  ${PYTHONHOME}/bin/python -v ./numpyTest.py
 EOF

 qsub numpyTest.pbs

Install ase/gpaw-setups (**Note**: use the latest releases)::

  cd ${sw_home}
  wget --no-check-certificate https://wiki.fysik.dtu.dk/ase-files/python-ase-3.1.0.846.tar.gz
  tar zxf python-ase-3.1.0.846.tar.gz
  wget --no-check-certificate http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.5.3574.tar.gz
  tar zxf gpaw-setups-0.5.3574.tar.gz

  cp -r python-ase-3.1.0.846 gpaw-setups-0.5.3574 /tmp/work/$USER
  cd /tmp/work/$USER
  ln -s python-ase-3.1.0.846 ase

Install gpaw (**Note**: instructions valid from the **5232** release)::

  cd ${sw_home}
  wget --no-check-certificate https://wiki.fysik.dtu.dk/gpaw/gpaw-0.7.5232.tar.gz
  tar zxf gpaw-0.7.5232.tar.gz
  cd gpaw-0.7.5232
  wget --no-check-certificate http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Cray/customize_jaguar.py -O customize.py
  ${PYTHON_DIR}/bin/python setup.py build_ext | tee build_ext.log
  cp -r ${sw_home}/gpaw-0.7.5232 /tmp/work/$USER
  cd /tmp/work/$USER
  ln -s gpaw-0.7.5232 gpaw

Test gpaw::

  cd /tmp/work/$USER

  cat <<EOF > ./gpawTest.pbs
  #!/bin/bash
  #PBS -l walltime=00:40:00,size=8
  #PBS -N gpawTest
  #PBS -A XXXXXX
  #PBS -j oe

  export PYTHONHOME=/tmp/work/$USER/Python-2.5.4-1
  export GPAW_SETUP_PATH=/tmp/work/$USER/gpaw-setups-0.5.3574
  export PYTHONPATH=/tmp/work/$USER/gpaw:/tmp/work/$USER/ase

  cd /tmp/work/$USER/gpaw/gpaw/test
  env | grep PYTHON
  env | grep LD_LIBRARY_PATH
  aprun -n4 /tmp/work/$USER/gpaw/build/bin.linux-x86_64-2.5/gpaw-python -v ./test.py
  EOF

  qsub gpawTest.pbs

