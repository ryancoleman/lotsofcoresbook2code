.. _Ubuntu_manual_install:

=============================
Manual installation on Ubuntu
=============================

Here you find information about the `Ubuntu distribution
<http://www.ubuntu.com/>`_.

Install ASE and the required packages, listed below, then
:ref:`download <download>` GPAW trunk or stable and modify
:file:`.bashrc` as detailed in the :ref:`installationguide`.  The
following instructions are for Ubuntu 8.10 or newer.  Instructions for
older versions can be found :ref:`further down the
page <ubuntu_oldversions>`.

**Warning** the versions of ScaLAPACK distributed by Ubuntu are old
and known to cause problems in GPAW.
See https://trac.fysik.dtu.dk/projects/gpaw/ticket/230 - use
ScaLAPACK at your own risk!

Required packages:

* python-dev
* python-numpy
* libopenblas-dev
* liblapack-dev
* libxc-dev

For Ubuntu/Debian versions where libxc is not available
install libxc manually, or download the Ubuntu
Trusty (14.04) deb packages and install with dpkg, e.g. for amd64::

   wget http://launchpadlibrarian.net/161692443/libxc-dev_2.0.2-1ubuntu1_amd64.deb
   wget http://launchpadlibrarian.net/161692442/libxc1_2.0.2-1ubuntu1_amd64.deb
   sudo dpkg -i libxc*_amd64.deb

Recommended:

* python-scientific
* python-matplotlib

Building documentation (XXX check whether this is precisely enough):

* python-sphinx
* texlive-latex-base
* texlive-latex-extra
* dvipng
* povray

To run parallel calculations, choose one of the three sections below.

Parallelization with OpenMPI
----------------------------

Install the packages

* openmpi-bin
* libopenmpi-dev

This will provide all parallel functionality except that of
BLACS/ScaLAPACK, which is normally not needed.

Parallelization with OpenMPI and BLACS/ScaLAPACK
------------------------------------------------

This is recommended for Ubuntu 9.10.  Install the package

* libscalapack-mpi-dev

Then build GPAW with the customize-file
:svn:`~doc/install/Linux/customize-ubuntu-sl-blacs-openmpi.py`:

.. literalinclude:: customize-ubuntu-sl-blacs-openmpi.py

Parallelization with MPI-LAM and BLACS/ScaLAPACK
------------------------------------------------

For Ubuntu 9.04 and possibly some older versions, there are no
BLACS/ScaLAPACK packages for use with OpenMPI.  In this case it is
recommended to use MPI-LAM.  Install the packages

* lam-runtime
* scalapack1-lam
* scalapack-lam-dev
* blacs1gf-lam
* blacsgf-lam-dev

Then build GPAW with the customize-file
:svn:`~doc/install/Linux/customize-ubuntu-sl-blacs-lam.py`:

.. literalinclude:: customize-ubuntu-sl-blacs-lam.py

.. _ubuntu_oldversions:

-----------------------
Ubuntu 8.04 or earlier
-----------------------

Install these packages:

* python-dev
* lapack3
* lapack3-dev
* refblas3
* refblas3-dev
* build-essential
* python-numpy
* python-numpy-ext

Recommended:

* atlas3-base
* atlas3-base-dev
* atlas3-headers
* python-scientific

GPAW will use atlas3 if available, which should increase performance.
Python-scientific is not strictly necessary, but some tests require
it.  Some packages in build-essential are likewise not necessary.
Note that this will only be a serial GPAW installation.
