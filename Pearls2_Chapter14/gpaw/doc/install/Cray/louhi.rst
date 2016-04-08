.. _louhi:

============================
louhi.csc.fi  (Cray XT4/XT5) 
============================

Here you find information about the the system
`<http://www.csc.fi/english/research/Computing_services/computing/servers/louhi>`_.

.. note::
   These instructions are up-to-date as of August 28th 2012.

GPAW
====

The recent operating system releases for Cray XT4/5 (CLE 2.2 UP01 and later) 
supports dynamic libraries which simplifies GPAW installation significantly.

These instructions for GPAW installation use Python 2.6.5 compiled
with GNU compiler suite, see the end of this page for instructions for
compiling Python.

First, load the Python module and set ``XTPE_LINK_TYPE`` environment
variable for dynamic linking::

  module load python
  module load hdf5-parallel
  setenv XTPE_LINK_TYPE dynamic

GPAW can now be build with a minimal ``customize.py``

.. literalinclude:: customize_louhi.py

Currently, there is a small bug in the Cray libraries which results in
failure when trying to build the ``_gpaw.so`` library. As the library
is not really needed in parallel calculations, the problem can be
circumvented by creating the file after running the setup.py script::

  python setup.py build_ext
  ...
  touch build/lib.linux-x86_64-2.6/_gpaw.so
  python setup.py install --home=...

Python and Numpy
================

Python can be compiled with PGI compiler as follows::

  setenv XTPE_LINK_TYPE dynamic
  ./configure --prefix=path_to_install CC=cc CXX=cc OPT=-fastsse LINKFORSHARED=-Wl,--export-dynamic
  make install

In order to use optimized BLAS with Numpy one has to first build a
CBLAS which is linked with Cray's optimized BLAS routines. First,
download the CBLAS source from netlib::

    wget http://www.netlib.org/blas/blast-forum/cblas.tgz
    tar -xzf cblas.tgz

Change to the CBLAS directory and copy ``Makefile.LINUX`` to
``Makefile.in``. Add correct compiler commands and paths to
``Makefile.in``::

  ...
  PLAT = louhi

  #-----------------------------------------------------------------------------
  # Libraries and includs
  #-----------------------------------------------------------------------------

  BLLIB =
  CBDIR = $(HOME)/CBLAS
  CBLIBDIR = $(CBDIR)/lib/$(PLAT)
  CBLIB = $(CBLIBDIR)/libcblas.a

  #-----------------------------------------------------------------------------
  # Compilers
  #-----------------------------------------------------------------------------

  CC = cc
  FC = ftn

  LOADER = $(FC)

  #-----------------------------------------------------------------------------
  # Flags for Compilers
  #-----------------------------------------------------------------------------

  CFLAGS = -O3 -DADD_ -fPIC
  FFLAGS = -O3 -fPIC

  ...

Finally, build CBLAS::

  make alllib

You are now ready to build Numpy with the newly created CBLAS
library. The standard Numpy tries to use only the ATLAS BLAS, and in
order to use different BLAS one has to manually edit the file
``numpy/core/setup.py``. Comment out an if statement as follows::

      def get_dotblas_sources(ext, build_dir):
        if blas_info:
            # if ('NO_ATLAS_INFO',1) in blas_info.get('define_macros',[]):
            #     return None # dotblas needs ATLAS, Fortran compiled blas will not be sufficient.
            return ext.depends[:1]

Then, add the correct libraries and paths to the file ``site.cfg``::

  [blas]
  blas_libs = cblas
  library_dirs = /home/csc/jenkovaa/CBLAS/lib/louhi

  [lapack]
  lapack_libs = sci
  library_dirs = /opt/xt-libsci/10.3.8/pgi/lib

Now, one should be able to build Numpy as usual::

  python setup.py install




