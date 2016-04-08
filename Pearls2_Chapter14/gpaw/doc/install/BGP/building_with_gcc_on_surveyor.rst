.. _building_with_gcc_on_surveyor:

==========================
Building with gcc compiler
==========================

NumPy
=====

If you do not wish to build NumPy for yourself, you can use one of the following versions::

 /soft/apps/python/python-2.6-cnk-gcc/numpy-1.2.1/lib/python2.6/site-packages
 /soft/apps/python/python-2.6-cnk-gcc/numpy-1.3.0/lib/python2.6/site-packages

Choose your version of NumPy accordingly. NumPy 1.3.0 officially supports Python 2.6, NumPy 1.2.1
is available as a fall back for use with Python 2.6. **We highly recommend that you use the pre-built
NumPy 1.3.0 rather than building your own.**

The **0.3** version of gpaw uses Numeric `<https://svn.fysik.dtu.dk/projects/gpaw/tags/0.3/>`_.

Get the Numeric-24.2 (**only** if you want to run the **0.3** version of gpaw)
and do this::

  $ wget http://downloads.sourceforge.net/numpy/Numeric-24.2.tar.gz
  $ gunzip -c Numeric-24.2.tar.gz | tar xf -
  $ cd Numeric-24.2
  $ /bgsys/drivers/ppcfloor/gnu-linux/bin/python setup.py install --root=$HOME/Numeric-24.2-1

The latest version of gpaw uses numpy `<https://svn.fysik.dtu.dk/projects/gpaw/trunk/>`_.

To build an optimized numpy for the compute nodes, based on ``goto`` blas, save the :svn:`numpy-1.0.4-gnu.py.patch.powerpc-bgp-linux-gfortran`
patch file
(modifications required to get powerpc-bgp-linux-gfortran instead of
gfortran compiler),
the :svn:`numpy-1.0.4-system_info.py.patch.lapack_bgp_goto_esslbg` patch file (lapack
section configured to use ``lapack_bgp`` and
blas section to use ``goto``, ``cblas_bgp``, and ``esslbg``),
and the :svn:`numpy-1.0.4-site.cfg.lapack_bgp_goto_esslbg` file (contains paths to
``lapack_bgp``, ``goto``, ``esslbg`` , ``cblas_bgp``,
and xlf* related libraries).

**Note** that ``lapack_bgp`` and ``cblas_bgp`` are not available on
``surveyor/intrepid``, to build use instructions from
`<http://www.pdc.kth.se/systems_support/computers/bluegene/LAPACK-CBLAS/LAPACK-CBLAS-build>`_. Python
requires all libraries to have names like ``liblapack_bgp.a``, so
please make the required links for ``lapack_bgp.a``, and
``cblas_bgp.a``. Moreover numpy requires that ``lapack_bgp``,
``goto``, ``esslbg``, and ``cblas_bgp`` reside in the same directory,
so choose a directory and edit
``numpy-1.0.4-site.cfg.lapack_bgp_goto_esslbg`` to reflect your
installation path (in this example
/home/dulak/from_Nils_Smeds/CBLAS_goto/lib/bgp). Include the directory
containing cblas.h in include_dirs. Change the locations of the
libraries to be used in the makefiles: /soft/apps/LIBGOTO and
/opt/ibmcmp/lib/bg.

**Warning** : If NumPy built using these libraries fails
with errors of kind "R_PPC_REL24 relocation at 0xa3d664fc for symbol sqrt"
- please add ``-qpic`` to compile options for both ``lapack_bgp`` and ``cblas_bgp``. 
After building ``lapack_bgp`` and ``cblas_bgp``, get numpy-1.0.4 and do this::

  $ wget http://downloads.sourceforge.net/numpy/numpy-1.0.4.tar.gz
  $ gunzip -c numpy-1.0.4.tar.gz | tar xf -
  $ mv numpy-1.0.4 numpy-1.0.4.optimized; cd numpy-1.0.4.optimized
  $ patch -p1 < ../numpy-1.0.4-gnu.py.patch.powerpc-bgp-linux-gfortran
  $ patch -p1 < ../numpy-1.0.4-system_info.py.patch.lapack_bgp_goto_esslbg
  $ cp ../numpy-1.0.4-site.cfg.lapack_bgp_goto_esslbg site.cfg
  $ ldpath=/bgsys/drivers/ppcfloor/gnu-linux/lib
  $ ldflags="-Wl,--allow-multiple-definition -L/opt/ibmcmp/xlmass/bg/4.4/bglib"
  $ root=$HOME/numpy-1.0.4-1.optimized
  $ p=/bgsys/drivers/ppcfloor/gnu-linux/bin/python
  $ c="\"/bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gcc -DNO_APPEND_FORTRAN -L/opt/ibmcmp/xlmass/bg/4.4/bglib\""
  $ MATHLIB="mass" LDFLAGS="$ldflags" LD_LIBRARY_PATH="$ldpath" CC="$c" $p setup.py install --root="$root"

NumPy built in this way does contain the
:file:`$root/bgsys/drivers/ppcfloor/gnu-linux/lib/python2.5/site-packages/numpy/core/_dotblas.so`
, and running the following python script results
in better time than the standard version of numpy (~156 vs. ~329 sec)
for ``numpy.dot`` operation (:svn:`~doc/install/BGP/numpy_dot.py`):

.. literalinclude:: numpy_dot.py

To build standard numpy, save the :svn:`~doc/install/BGP/numpy-1.0.4-gnu.py.patch` patch file
(modifications required to get mpif77 instead of gfortran compiler),
get and numpy-1.0.4 and do this::

  $ wget http://downloads.sourceforge.net/numpy/numpy-1.0.4.tar.gz
  $ gunzip -c numpy-1.0.4.tar.gz | tar xf -
  $ cd numpy-1.0.4
  $ patch -p1 < ../numpy-1.0.4-gnu.py.patch
  $ ldpath=/bgsys/drivers/ppcfloor/gnu-linux/lib
  $ root=$HOME/numpy-1.0.4-1
  $ p=/bgsys/drivers/ppcfloor/gnu-linux/bin/python
  $ c="\"mpicc\""
  $ LD_LIBRARY_PATH="$ldpath" CC="$c" $p setup.py install --root="$root"

The instructions for standard numpy also work for NumPy 1.2.1 and NumPy 1.3.0.

Build python-nose::

  $ wget http://python-nose.googlecode.com/files/nose-0.11.0.tar.gz
  $ tar zxf nose-0.11.0.tar.gz
  $ cd nose-0.11.0
  $ p=/bgsys/drivers/ppcfloor/gnu-linux/bin/python
  $ $p setup.py install --root=${HOME}/python-nose-0.11.0-1 2>&1 | tee install.log

GPAW
====

Step 1
======

Download all the necessary packages:

- `ase <https://wiki.fysik.dtu.dk/ase/download.html#latest-stable-release>`_

- `gpaw <https://wiki.fysik.dtu.dk/gpaw/download.html#latest-stable-release>`_

- `gpaw-setups <https://wiki.fysik.dtu.dk/gpaw/setups/setups.html>`_

Step 2
======

Add the following environment variables

.. literalinclude:: surveyor.softenvrc

to your own :file:`.softenvrc` and type::

  resoft

to update your environment in the main login terminal.

Step 3
======

This step is obsolete with revision 6116.

Because the ``popen3`` function is missing, you will need to remove all the
contents of the :file:`gpaw/version.py` file
after ``ase_required_svnversion =``.
The same holds for :file:`ase/version.py` in the ase installation!
Suggestions how to skip the ``popen3`` testing in
:file:`gpaw/version.py` on BG/P are welcome!

Step 4
======

A number of the GPAW source files in ``gpaw/c`` directory are built using
the ``distutils`` module which makes it difficult to control the flags
which are passed to the gnu compiler. A workaround is to use the following python
script: :svn:`bgp_gcc.py`.  Additionaly, it is
desirable to static link as many libraries as possible. This requires
bypassing the mpi wrapper to the compiler using another python
script :svn:`bgp_gcc_linker.py`.
Lastly, we must use these work arounds in conjuction with two
configures files :svn:`customize_surveyor_gcc.py`  and
:svn:`config_surveyor.py`, the latter requires
renaming to ``config.py`` in the top level directory.

.. literalinclude:: bgp_gcc.py

.. literalinclude:: bgp_gcc_linker.py

.. literalinclude:: customize_surveyor_gcc.py

Download these scripts into the top level GPAW directory::

  export GPAW_TRUNK=http://svn.fysik.dtu.dk/projects/gpaw/trunk
  wget --no-check-certificate ${GPAW_TRUNK}/doc/install/BGP/bgp_gcc.py
  chmod u+x bgp_gcc.py
  wget --no-check-certificate ${GPAW_TRUNK}/doc/install/BGP/bgp_gcc_linker.py
  chmod u+x bgp_gcc_linker.py
  wget --no-check-certificate ${GPAW_TRUNK}/doc/install/BGP/customize_surveyor_gcc.py
  wget --no-check-certificate ${GPAW_TRUNK}/doc/install/BGP/config_surveyor.py
  mv config_surveyor.py config.py

Finally, we build GPAW by typing::

  /bgsys/drivers/ppcfloor/gnu-linux/bin/python setup.py build_ext --ignore-numpy --customize=customize_surveyor_gcc.py 2>&1 | tee build_ext.log

If an optimized version of NumPy is in your $PYTHONPATH you may need append  "--ignore-numpy".

Additional BG/P specific hacks
===============================

A FLOPS (floating point per second) counter and a number of other hardware counters can be enabled with the macro::

  define_macros += [('GPAW_HPM',1)]

This hpm library is available on the BG/P machines at Argonne National Laboratory. It will produce two files for each core: ``hpm_flops.$rank`` and ``hpm_data.$rank``. The latter one contains a number of additional hardware counters. There are four cores per chip and data for only two of the four cores can be collected simultaneously. This is set through an environment variable which is passed to Cobalt with the *--env*  flag. *BGP_COUNTER_MODE=0* specifies core 1 and 2, while *BGP_COUNTER_MODE=1* specifies core 3 and 4. 

A mapfile for the ranks can be generated by adding another macro to customize.py::
  
  define_macros += [('GPAW_MAP',1)]


Submitting jobs
==================

This is an example submission script
:svn:`surveyor.sh` for use with the Cobalt scheduler:

.. literalinclude:: surveyor.sh

Users should read the :ref:`parallel_runs` page and the 
:ref:`bgp_performance` page.

Running from ramdisk
======================

There are two I/O bottlenecks present when running GPAW at scale: 1)
loading .py(c) and 2) loading .so. For jobs requiring 8192 vn nodes
or larger, the initialization time was measured at 40+ minutes.

A work-around developed for the BlueGene/P system at Argonne National 
Laboratory was to install GPAW and all supporting libraries on a
ramdisk. The supporting libraries include::

   Python 2.6 
   NumPy 1.3.0
   ASE
   GPAW
   <other system shared libraries>

The ramdisk utility provided IBM has some weird quirks that one 
should be aware of: 1) files cannot be installed in the root
of the ramdisk, everything should go in a directory like ``/gpaw``
2) empty directories and their child directores are not recognized
properly;  the solution is to create a zero-size file, e.g. ``touch
README``, in each empty directory. 3) symbolic links are not supported

The top-level directory/file structure should look something like this::

  $HOME/bgpramdisk/Makefile <file>
  $HOME/bgpramdisk/fs
  $HOME/bgpramdisk/fs/gpaw
  $HOME/bgpramdisk/fs/gpaw/README <file>
  $HOME/bgpramdisk/fs/gpaw/lib
  $HOME/bgpramdisk/fs/gpaw/V1R4M2

An example Makefile is provided :svn:`Makefile.ramdisk`, note that
BGP_LINUX_OS_PATH version will be a function of the OS driver. The
group and installation directory will be provided by the BlueGene/P sys admin.

Copy all the supporting libraries unto the ramdisk. For
simplicity, we install the Python packages into the site-packages
directory which is automatically searched by the Python Intrepreter::

  cp -r /bgsys/driver/ppcfloor/gnu/linux/lib/python2.6 $HOME/bgpramdisk/fs/gpaw/lib/
  cp -r /soft/apps/python/python-2.6-cnk-gcc/numpy-1.3.0/lib/site-packages/numpy $HOME/bgpramdisk/fs/gpaw/lib/python2.6/site-packages/
  cp -r $HOME/gpaw-<version>/gpaw  $HOME/bgpramdisk/fs/gpaw/lib/python2.6/site-packages/
  cp -r $HOME/ase-<version>/ase  $HOME/bgpramdisk/fs/gpaw/python2.6/lib/site-packages/

Make sure that all these libraries are byte compiled. Python 2.6 and
NumPy 1.3.0 should already be byte compiled, but it is likely that ASE
and GPAW won't be. ``cd`` into their respective directories on the
ramdisk and type::

   /bgsys/drivers/ppcfloor/gnu-linux/bin/python -m compileall .

In order to save space on the ramdisk, delete the ``*.py``,  but keep the ``*.pyc.``::
  
   find -name "*.py" -exec rm -rf {} \;

The MPI shared libraries can also installed on the ramdisk so they will
seamlessly work with the TAU Performance System, but static libraries 
can be used instead::
  
  cp /bgsys/drivers/ppcfloor/comm/default/lib/libmpich.cnk.so.1.1 $HOME/bgpramdisk/fs/gpaw/V1R4M2/
  cp /bgsys/drivers/ppcfloor/comm/sys/lib/libdcmfcoll.cnk.so $HOME/bgpramdisk/fs/gpaw/V1R4M2/
  cp /bgsys/drivers/ppcfloor/comm/sys/lib/libdcmf.cnk.so $HOME/bgpramdisk/fs/gpaw/V1R4M2/
  cp /bgsys/drivers/ppcfloor/runtime/SPI/libSPI.cna.so  $HOME/bgpramdisk/fs/gpaw/V1R4M2/

TAU profiling also requires some additional files for automatic and
manual instrumentation::

  cp /soft/apps/tau/tau-2.19.2/bgp/lib/ltau.py $HOME/bgpramdisk/fs/gpaw/tau-2.19.2/
  cp /soft/apps/tau/tau-2.19.2/bgp/lib/tau.py $HOME/bgpramdisk/fs/gpaw/tau-2.19.2/  
  cp /soft/apps/tau/tau-2.19.2/bgp/lib/ibTAUsh-phase-bgptimers-gnu-mpi-python-pdt.so $HOME/bgpramdisk/fs/gpaw/tau-2.19.2/
  cp /soft/apps/tau/tau-2.19.2/bgp/lib/ibTAUsh-phase-bgptimers-gnu-mpi-python-pdt.so $HOME/bgpramdisk/fs/gpaw/tau-2.19.2/ctau_impl.so
  cp /soft/apps/tau/tau-2.19.2/bgp/lib/ibTAUsh-phase-bgptimers-gnu-mpi-python-pdt.so $HOME/bgpramdisk/fs/gpaw/tau-2.19.2/pytau.so

Lastly, there will be some changes to the environment variables in your submission script::

  PYTHONHOME=/gpaw
  LD_LIBRARY_PATH=/lib:/gpaw/V1R4M2

PYTHONPATH should be empty unless you have installed another software
package on the ramdisk. In any case, it should not point to any physical diskspace
