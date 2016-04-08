.. _todolist:

=========
Todo list
=========

.. seealso::

   * GPAW's `issue tracker <http://trac.fysik.dtu.dk/projects/gpaw/report/1>`__
   * ASE's `issue tracker <http://trac.fysik.dtu.dk/projects/ase/report/1>`__


Simple tasks
============

* Remove this nonsense from GPAW's text output: "Bands to Converge:
  -20 Lowest Bands".

* Merge ``Domain`` and ``GridDescriptor`` classes into one class and
  remove ``xxxucell_cv`` and ``xxxiucell_cv`` members.

* Remove ``get_scaled_positions() % 1`` statements.  Should not be
  needed since ASE version 3.4.

* The :ref:`pawxml` needs to be extended to include MGGA quantities
  like the core kinetic energy density.

* Write dipole moment to trajectory files.

* "Broken symmetry": We need a better error message.  This happens
  often when ASE's vibration module tries to move atoms.

* Calculate a memory estimate for vdW-DF calculations.

* Stop using ``wfs.ibzk_qc`` and other things that are in ``wfs.kd``.

* Use ``ase.Atoms.get_charges`` instead of the ``charge`` keyword.


Questions
=========

* ``DeprecationWarning``'s have been silenced in Python 2.7.  What do we do?

* Should we warn users if they use a ridiculously slow BLAS library?

* Can we use ctypes to link to MPI, libxc, and GPAW's own C-code?
  Would this work on BlueGene and Cray?


Larger projects
===============

* Allow the number of grid points to be multiples of 3, 5, 7, 11 and so on.

* Rename gpaw.utilities to gpaw.utils as in ASE and clean up!

* Invent more flexible specification of setups using file names or
  directories: ``setups={'La': '~/mysetups'}``.

* Finish dipole layer implementation.

* Update our copy of the libxc source code and integrate our MGGA changes.

* Implement the HSE03 hybrid functional.


Documentation
=============

* Write a Sphinx-role that can convert a `DOI number
  <http://dx.doi.org>`_ to title, authers and journal.

* Remove png's from svn-repository.

* Update :ref:`The big picture <the_big_picture>`: remove the two
  ``gd`` boxes and the ``kpt_comm`` box.

* Write a C-style guide.


Important setup problems to fix
===============================

* Improve Ru setup.  Problem with nitrogen adsorption energy on
  Ru(0001) surface: Improved f-scattering, do we need 4p semi-core
  states? This problem has been fixed in the 0.9 release of the
  setups and covered by the following test :svn:`~gpaw/test/big/Ru001`.
* Improve Mn setup.  Problem with strange states in the band-gap for
  AFM-II MnO.
* Fine-tune Fourier filtering of projector functions.  There is still
  room for improvements.
