.. _rbgc:

=======================
bcssh.rochester.ibm.com
=======================

Instructions below are valid for ``frontend-13`` and the filesystem
:file:`/gpfs/fs2/frontend-13`.

The latest version of gpaw uses numpy
`<https://svn.fysik.dtu.dk/projects/gpaw/trunk/>`_.

To build an optimized (consider to build based on ``goto`` blas to achieve the best performance
see :ref:`building_with_gcc_on_surveyor`) numpy,
save the :svn:`~doc/install/BGP/numpy-1.0.4-gnu.py.patch.powerpc-bgp-linux-gfortran` patch file
(modifications required to get powerpc-bgp-linux-gfortran instead of
gfortran compiler),
the :svn:`~doc/install/BGP/numpy-1.0.4-system_info.py.patch.lapack_bgp_esslbg` patch file (lapack
section configured to use ``lapack_bgp`` and
blas section to use ``esslbg`` and ``cblas_bgp``),
and the :svn:`~doc/install/BGP/numpy-1.0.4-site.cfg.lapack_bgp_esslbg` file (contains paths to
``lapack_bgp``, ``esslbg`` , ``cblas_bgp``, and xlf* related libraries).

**Note** that ``lapack_bgp`` and ``cblas_bgp`` are not available on
``frontend-13``, to build use instructions from
`<http://www.pdc.kth.se/systems_support/computers/bluegene/LAPACK-CBLAS/LAPACK-CBLAS-build>`_. Python
requires all librairies to have names like ``liblapack_bgp.a``, so
please make the required links for ``lapack_bgp.a`` and
``cblas_bgp.a``. Moreover numpy requires that ``lapack_bgp``,
``esslbg``, and ``cblas_bgp`` reside in the same directory, so choose
a directory and edit ``numpy-1.0.4-site.cfg.lapack_bgp_esslbg`` to
reflect your installation path (in this example
/home/dulak/from_Nils_Smeds/CBLAS/lib/bgp). Include the directory
containing cblas.h in include_dirs. These instructions are valid
also for Surveyor/Intrepid with the following locations of the
libraries to be used in the makefiles: /soft/apps/ESSL-4.4/lib and
/opt/ibmcmp/lib/bg.

**Warning** - if numpy built using these libraries fails
with errors of kind "R_PPC_REL24 relocation at 0xa3d664fc for symbol sqrt"
- please add ``-qpic`` to compile options for both ``lapack_bgp`` and ``cblas_bgp``. 
After bulding ``lapack_bgp`` and ``cblas_bgp``, get numpy-1.0.4 and do this::

  $ wget http://downloads.sourceforge.net/numpy/numpy-1.0.4.tar.gz
  $ gunzip -c numpy-1.0.4.tar.gz | tar xf -
  $ mv numpy-1.0.4 numpy-1.0.4.optimized; cd numpy-1.0.4.optimized
  $ patch -p1 < ../numpy-1.0.4-gnu.py.patch.powerpc-bgp-linux-gfortran
  $ patch -p1 < ../numpy-1.0.4-system_info.py.patch.lapack_bgp_esslbg
  $ cp ../numpy-1.0.4-site.cfg.lapack_bgp_esslbg site.cfg
  $ ldpath=/bgsys/drivers/ppcfloor/gnu-linux/lib
  $ mkdir /gpfs/fs2/frontend-13/$USER
  $ root=/gpfs/fs2/frontend-13/$USER/numpy-1.0.4-1.optimized
  $ p=/bgsys/drivers/ppcfloor/gnu-linux/bin/python
  $ c="\"/bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gcc -DNO_APPEND_FORTRAN\""
  $ LD_LIBRARY_PATH="$ldpath" CC="$c" $p setup.py install --root="$root"

Numpy built in this way does contain the
:file:`$root/bgsys/drivers/ppcfloor/gnu-linux/lib/python2.5/site-packages/numpy/core/_dotblas.so`
, but running the following python
script (save it as :file:`/gpfs/fs2/frontend-13/$USER/numpy_dot.py`) results
in the same time as for the standard version of numpy (~329 sec)
for ``numpy.dot`` operation (:svn:`~doc/install/BGP/numpy_dot.py`):

.. literalinclude:: numpy_dot.py

Use the following command to submit this job ``cd
/gpfs/fs2/frontend-13/$USER; llsubmit numpy.llrun``, with the
following :svn:`~doc/install/BGP/numpy.llrun` file:

.. literalinclude:: numpy.llrun

**Note** the colon before and after the string when setting pythonpath!

Here is how you build the standard numpy::

  $ gunzip -c numpy-1.0.4.tar.gz | tar xf -
  $ cd numpy-1.0.4
  $ patch -p1 < ../numpy-1.0.4-gnu.py.patch.powerpc-bgp-linux-gfortran
  $ ldpath=/bgsys/drivers/ppcfloor/gnu-linux/lib
  $ mkdir /gpfs/fs2/frontend-13/$USER
  $ root=/gpfs/fs2/frontend-13/$USER/numpy-1.0.4-1
  $ p=/bgsys/drivers/ppcfloor/gnu-linux/bin/python
  $ c="\"/bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gcc\""
  $ LD_LIBRARY_PATH="$ldpath" CC="$c" $p setup.py install --root="$root"

Build GPAW
(``PYTHONPATH=/gpfs/fs2/frontend-13/mdulak/numpy-1.0.4-1.optimized/bgsys/drivers/ppcfloor/gnu-linux/lib/python2.5/site-packages
LD_LIBRARY_PATH="$ldpath" $p setup.py build_ext --do-not-force-inclusion-of-numpy``) in
:file:`/gpfs/fs2/frontend-13/$USER/gpaw` (you need to install the ase
also somewhere below :file:`/gpfs/fs2/frontend-13/$USER`!)  with this
:svn:`~doc/install/BGP/customize_rbgc.py` file:

.. literalinclude:: customize_rbgc.py

Because of missing ``popen3`` function you need to remove all the
contents of the :file:`gpaw/version.py` file after ``version =
'0.4'``.  The same holds for :file:`ase/version.py` in the ase
installation!  Suggestions how to skip the ``popen3`` testing in
:file:`gpaw/version.py` on BGP are welcome!

Note that only files located below :file:`/gpfs/fs2/frontend-13` are
accesible to the compute nodes (even python scripts!).  A gpaw script
:file:`/gpfs/fs2/frontend-13/$USER/gpaw/test/CH4.py` can be submitted to
32 CPUs in the single mode (SMP) for 30 minutes using `LoadLeveler
<http://www.fz-juelich.de/jsc/ibm-bgl/usage/loadl/>`_ like this::

  cd /gpfs/fs2/frontend-13/$USER
  llsubmit gpaw-script.llrun

where :svn:`~doc/install/BGP/gpaw-script.llrun` looks like this:

.. literalinclude:: gpaw-script.llrun

Absolute paths are important!

It's convenient to customize as in :file:`gpaw-qsub.py` which can be
found at the :ref:`parallel_runs` page.
