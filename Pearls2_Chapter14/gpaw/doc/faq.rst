.. _faq:

==========================
Frequently Asked Questions
==========================

General
=======

Citation: how should I cite GPAW?
---------------------------------

If you find GPAW useful in your research please cite the original reference:

   | J. J. Mortensen, L. B. Hansen , and K. W. Jacobsen
   | `Real-space grid implementation of the projector augmented wave method`__
   | Physical Review B, Vol. **71**, 035109, 2005
  
   __ http://dx.doi.org/10.1103/PhysRevB.71.035109

and the major GPAW review:

   | J. Enkovaara, C. Rostgaard, J. J. Mortensen et al.
   | `Electronic structure calculations with GPAW: a real-space implementation of the projector augmented-wave method`__ 
   | J. Phys.: Condens. Matter **22**, 253202 (2010)

   __ http://stacks.iop.org/0953-8984/22/253202


together with :ase:`ASE <>` citation
(see :ase:`Citation: how should I cite ASE?<faq.html>`).

If you are using the time-dependent DFT part of the code, please cite also:

   | M. Walter, H. Häkkinen, L. Lehtovaara, M. Puska, J. Enkovaara, C. Rostgaard and J. J. Mortensen
   | `Time-dependent density-functional theory in the projector augmented-wave method`__
   | Journal of Chemical Physics, Vol. **128**, 244101, 2008

   __ http://link.aip.org/link/?JCP/128/244101

If you use the :ref:`localized basis set <lcao>`, please cite also:

   | A. H. Larsen, M. Vanin, J. J. Mortensen, K. S. Thygesen, and K. W. Jacobsen
   | `Localized atomic basis set in the projector augmented wave method`__
   | Physical Review B, Vol. **80**, 195112, 2009

   __ http://dx.doi.org/10.1103/PhysRevB.80.195112
   
If you use the :ref:`df_tutorial`, please cite also:

   | Jun Yan, Jens. J. Mortensen, Karsten W. Jacobsen, and Kristian S. Thygesen
   | `Linear density response function in the projector augmented wave method: Applications to solids, surfaces, and interfaces`__
   | Physical Review B Vol. **83**, 245122, 2011

   __ http://link.aps.org/doi/10.1103/PhysRevB.83.245122

If you use the :ref:`gw_tutorial`, please cite also:

   | F. Hüser, T. Olsen, and K. S. Thygesen
   | `Quasiparticle GW calculations for solids, molecules, and two-dimensional materials`__
   | Physical Review B Vol. **87**, 235132, 2013

   __ http://link.aps.org/doi/10.1103/PhysRevB.87.235132

BibTex (:svn:`doc/GPAW.bib`):

.. literalinclude:: GPAW.bib


How do you pronounce GPAW?
--------------------------

In English: "geepaw" with a long "a".

In Danish: Først bogstavet "g", derefter "pav": "g-pav".

In Finnish: supisuomalaisittain "kee-pav".

In Polish: "gyeh" jak `"Gie"rek <http://en.wikipedia.org/wiki/Edward_Gierek>`_, "pav" jak `paw <http://pl.wikipedia.org/wiki/Paw_indyjski>`_: "gyeh-pav".

Download
========

Trying to checkout the code via SVN resulted::

 [~]$ svn checkout "https://svn.fysik.dtu.dk/projects/gpaw/trunk"
 svn: Unrecognized URL scheme 'https://svn.fysik.dtu.dk/projects/gpaw/trunk'

This error is diplayed in case the library 'libsvn_ra_dav' is missing
on your system. The library is used by SVN, but is not installed by
default.


Compiling the C-code
====================

For architecture dependent settings see the
:ref:`platforms_and_architectures` page.

Compilation of the C part failed::

 [~]$ python2.4 setup.py build_ext
 building '_gpaw' extension
 pgcc -fno-strict-aliasing -DNDEBUG -O2 -g -pipe -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -m64 -D_GNU_SOURCE -fPIC -fPIC -I/usr/include/python2.4 -c c/localized_functions.c -o build/temp.linux-x86_64-2.4/c/localized_functions.o -Wall -std=c99
 pgcc-Warning-Unknown switch: -fno-strict-aliasing
 PGC-S-0040-Illegal use of symbol, _Complex (/usr/include/bits/cmathcalls.h: 54)

You are probably using another compiler, than it was used for
compiling python. Undefine the environment variables CC, CFLAGS and
LDFLAGS with::

 # sh/bash users:
 unset CC; unset CFLAGS; unset LDFLAGS
 # csh/tcsh users: 
 unsetenv CC; unsetenv CFLAGS; unsetenv LDFLAGS

and try again.

Calculation does not converge
=============================

Consult the :ref:`convergence` page.

Poisson solver did not converge!
================================

If you are doing a spin-polarized calculation for an isolated molecule, 
then you should set the Fermi temperature to a low value.

You can also try to set the number of grid points to be divisible by 8. 
Consult the :ref:`poisson_performance` page.

How to switch between several GPAW versions
===========================================

For each GPAW installation use a separate, modified submit tool:
:svn:`~doc/documentation/parallel_runs/gpaw-qsub`.

Assuming that your :ref:`developer_installation` is under
:file:`~/gpaw.test`, and the :command:`gpaw-python` under
:file:`~/gpaw.test/build/bin.linux-x86_64-2.3/`, modify the submit
tool: :svn:`~doc/documentation/parallel_runs/gpaw-qsub`:

* set the :envvar:`PYTHONPATH` and :envvar:`PATH` passed to :command:`mpirun`::

   ...
   'export PYTHONPATH=${HOME}/gpaw.test:${PYTHONPATH} && ' +
   'export PATH=${HOME}/gpaw.test/build/bin.linux-x86_64-2.3:${PATH} && ' +
   'mpirun')

* make sure that the corresponding :command:`gpaw-python` is used::

   os.system('%s gpaw-python JOB' % (mpirun))

Alternatively, instead of modifying
:svn:`~doc/documentation/parallel_runs/gpaw-qsub`
create a bash function - see :ref:`Niflheim` for details.

Tests fail!
===========

Please report the failing test as described on :ref:`running_tests`.

