.. _SUNCAT:

======
SUNCAT
======

(NOTE: With MKL 10.3 we have seen hangs in the early mpireduce calls for
a small number of calculations.  Until I have understood this I am
backing out to MKL 10.2.)

At SLAC we compiled GPAW for RHEL5 x86_64, on intel Xeon 5650 with intel compilers and mkl. This improved the 8-core performance benchmark by 13% compared to the opencc/ACML approach.

================ ==================
Package          Version
================ ==================
python           2.4
gpaw             0.8.0.7419 
ase              3.5.0.1919
numpy            1.4.1
openmpi          1.4.3
mkl              10.3
intel compilers  11.1 (includes mkl 10.2 by default)
================ ==================

openmpi
=======

openmpi was built with the intel compilers as follows::

   $ ./configure --prefix=/nfs/slac/g/suncatfs/sw/gpawv15/install CC=icc CXX=icpc F77=ifort FC=ifort
   $ make
   $ make install

numpy
=====

Build in usual fashion. At the moment we use default gnu compilers for numpy, since gpaw performance benchmark drops by 3% when it is built with icc/mkl/dotblas, for reasons that are not understood. Also, some gpaw self-tests start to fail.

gpaw
====

For this we use
:svn:`~doc/install/Linux/SUNCAT/customize_mkl10.3.py`:

.. literalinclude:: customize_mkl10.3.py

Note that this customize.py works only with MKL version 10.3 which has simplified linking.

The environment settings (valid at SUNCAT) to be able to link and run:

.. literalinclude:: setupenv

MKL 10.2 Notes
==============

For historical reasons, we also include the customize.py for MKL 10.2:

.. literalinclude:: customize_mkl10.2.py

This older version requires a fairly bad hack to make it work in all cases::

   $ setenv LD_PRELOAD libmkl_core.so:libmkl_sequential.so

I believe this is because python uses "dlopen" for shared libraries, which has troubles with the circular dependencies present in MKL 10.2.

This hack can cause (ignorable) errors from unrelated commands like "ping" which prevents the use of LD_PRELOAD for security reasons.