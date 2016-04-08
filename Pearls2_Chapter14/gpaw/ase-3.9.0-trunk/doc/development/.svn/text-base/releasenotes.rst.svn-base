.. _releasenotes:

=============
Release notes
=============


Development version in trunk
============================

:trac:`trunk <>`.
          
* Genetic algorithm implemented; :mod:`ase.ga`. This can be used
  for the optimization of: atomic cluster structure, materials
  properties by use of template structures. Extension to other projects
  related to atomic simulations should be straightforward.

* The :func:`ase.lattice.bulk` function can now build the Wurtzite structure.

* The :class:`ase.utils.timer.Timer` was moved from GPAW to ASE.


Version 3.8.0
=============

22 October 2013: :trac:`tags/3.8.0 <../tags/3.8.0>`.

* ASE's :mod:`gui <gui>` renamed from ``ag`` to ``ase-gui``.
* New :ref:`STM <stm>` module.
* Python 2.6 is now a requirement.
* The old :func:`ase.structure.bulk` function is now deprecated.
  Use the new one instead (:func:`ase.lattice.bulk`).
* We're now using BuildBot for continous integration:
  https://ase-buildbot.fysik.dtu.dk/waterfall
* New interface to the JDFTx code.


Version 3.7.0
=============

13 May 2013: :trac:`tags/3.7.0 <../tags/3.7.0>`.

* ASE's GUI can now be configured to be more friendly to visually
  impaired users: :ref:`high contrast`.

* The :class:`ase.neb.NEB` object now accepts a list of spring constants.

* *Important backwards incompatible change*: The
  :func:`ase.lattice.surface.surface` function now returns a
  right-handed unit cell.

* Mopac, NWChem and Gaussian interfaces and EAM potential added.

* New :meth:`~ase.atoms.Atoms.set_initial_charges` and
  :meth:`~ase.atoms.Atoms.get_initial_charges` methods.  The
  :meth:`~ase.atoms.Atoms.get_charges` method will now ask the
  calculator to calculate the atomic charges.

* The :ref:`aep1` has been implemented and 6 ASE calculators are now
  based on the new base classes.

* ASE now runs on Windows and Mac.

* :ref:`mhtutorial` added to ASE.


Version 3.6.0
=============

24 Feb 2012: :trac:`tags/3.6.0 <../tags/3.6.0>`.

* ASE GUI translations added, available: da_DK, en_GB, es_ES.

* New function for making surfaces with arbitrary Miller indices with
  the smallest possible surface unit cell:
  ase.lattice.surface.surface()

* New ase.lattice.bulk() function.  Will replace old
  ase.structure.bulk() function.  The new one will produce a more
  natural hcp lattice and it will use experimental data for crystal
  structure and lattice constants if not provided explicitely.

* New values for ase.data.covalent_radii from Cordeo *et al.*.

* New command line tool: :ref:`command line tool` and tests based on it:
  abinit, elk, fleur, nwchem.

* New crystal builder for ase-gui

* Van der Waals radii in ase.data

* ASE's GUI (ase-gui) now supports velocities for both graphs and coloring

* Cleaned up some name-spaces:

  * ``ase`` now contains only :class:`~ase.atoms.Atoms` and
    :class:`~ase.atom.Atom`
  * ``ase.calculators`` is now empty


Version 3.5.1
=============

24 May 2011: :trac:`tags/3.5.1 <../tags/3.5.1>`.

* Problem with parallel vibration calculations fixed:
  `Ticket #80 <https://trac.fysik.dtu.dk/projects/ase/ticket/80>`_.


Version 3.5.0
=============

13 April 2011: :trac:`tags/3.5.0 <../tags/3.5.0>`.

* Improved EMT potential:  uses a
  :class:`~ase.calculators.neighborlist.NeighborList` object and is
  now ASAP_ compatible.

* :mod:`BFGSLineSearch <optimize.bfgslinesearch>` is now the default
  (``QuasiNewton==BFGSLineSearch``).

* There is a new interface to the LAMMPS molecular dynamics code.

* New :mod:`phonons` module.

* Van der Waals corrections for DFT, see GPAW_ usage.

* New :class:`~ase.io.bundletrajectory.BundleTrajectory` added.

* Updated GUI interface:

  * Stability and usability improvements.
  * Povray render facility.
  * Updated expert user mode.
  * Enabled customization of colours and atomic radii.
  * Enabled user default settings via :file:`~/.ase/gui.py`.

* :mod:`Database library <data>` expanded to include:
  
  * The s22, s26 and s22x5 sets of van der Waals bonded dimers and
    complexes by the Hobza group.
  * The DBH24 set of gas-phase reaction barrier heights by the Truhlar
    group.

* Implementation of the Dimer method.


.. _ASAP: http://wiki.fysik.dtu.dk/asap
.. _GPAW: https://wiki.fysik.dtu.dk/gpaw/documentation/xc/vdwcorrection.html


Version 3.4.1
=============

11 August 2010: :trac:`tags/3.4.1 <../tags/3.4.1>`.
