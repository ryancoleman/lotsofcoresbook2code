.. _releasenotes:

=============
Release notes
=============


Development version in trunk
============================

:trac:`trunk <>`.

* New :ref:`symmetry <manual_symmetry>` keyword.  Replaces ``usesymm``.

* Use non-symmorphic symmetries: combining fractional translations with
  rotations, reflections and inversion.  Use
  ``symmetry={'symmorphic': False}`` to turn this feature on.

* New :ref:`forces <manual_convergence>` keyword in convergence.  Can
  be used to calculate forces to a given precision.

Version 0.10.0
==============

8 April 2014: :trac:`tags/0.10.0 <../tags/0.10.0>`.

* Corresponding ASE release: ase-3.8.1_

* Default eigensolver is now the Davidson solver.

* Default density mixer parameters have been changed for calculations
  with periodic boundary conditions.  Parameters for that case:
  ``Mixer(0.05, 5, 50)`` (or ``MixerSum(0.05, 5, 50)`` for spin-paired
  calculations.  Old parameters: ``0.1, 3, 50``.
  
* Default is now ``occupations=FermiDirac(0.1)`` if a
  calculation is periodic in at least one direction,
  and ``FermiDirac(0.0)`` otherwise (before it was 0.1 eV for anything
  with **k**-points, and 0 otherwise).

* Calculations with a plane-wave basis set are now officially supported.

* :ref:`One-shot GW calculations <gw_theory>` with full frequency
  integration or plasmon-pole approximation.
  
* Beyond RPA-correlation: `using renormalized LDA and PBE
  <https://trac.fysik.dtu.dk/projects/gpaw/browser/branches/sprint2013/doc/tutorials/fxc_correlation>`_.

* :ref:`bse`.

* Improved RMM-DIIS eigensolver.

* Support for new libxc 2.0.1.  libxc must now be built separately from GPAW.

* MGGA calculations can be done in plane-wave mode.

* Calculation of the stress tensor has been implemented for plane-wave
  based calculation (except MGGA).

* MGGA: number of neighbor grid points to use for FD stencil for
  wave function gradient changed from 1 to 3.

* New setups: Y, Sb, Xe, Hf, Re, Hg, Tl, Rn

* Non self-consistent calculations with screened hybrid functionals
  (HSE03 and HSE06) can be done in plane-wave mode.

* Modified setups:

  .. note::

     Most of the new semicore setups currently require
     :ref:`eigensolver <manual_eigensolver>` ``dav``, ``cg``
     eigensolvers or ``rmm-diis`` eigensolver with a couple of iterations.

  - improved eggbox: N, O, K, S, Ca, Sc, Zn, Sr, Zr, Cd, In, Sn, Pb, Bi

  - semicore states included: Na, Mg, V, Mn, Ni,
    Nb, Mo, Ru (seems to solve the Ru problem :trac:`gpaw/test/big/Ru001`),
    Rh, Pd, Ag, Ta, W, Os, Ir, Pt

  - semicore states removed: Te

  - elements removed: La (energetics was wrong: errors ~1eV per unit cell
    for PBE formation energy of La2O3 wrt. PBE benchmark results)

  .. note::

     For some of the setups one has now a choice of different
     number of valence electrons, e.g.::

       setups = {'Ag': '11'}

     See :ref:`manual_setups` and list the contents of :envvar:`GPAW_SETUP_PATH`
     for available setups.

* new ``dzp`` basis set generated for all the new setups, see
  https://trac.fysik.dtu.dk/projects/gpaw/ticket/241


Version 0.9.0
=============

7 March 2012: :trac:`tags/0.9.0 <../tags/0.9.0>`.

* Corresponding ASE release: ase-3.6_

* Convergence criteria for eigenstates changed: The missing volume per
  grid-point factor is now included and the units are now eV**2. The
  new default value is 4.0e-8 eV**2 which is equivalent to the old
  default for a grid spacing of 0.2 Å.

* GPAW should now work also with NumPy 1.6.

* Much improved :ref:`command line tool` now based on the `new
  tool`_ in ASE.


.. _new tool: https://wiki.fysik.dtu.dk/ase/ase/cmdline.html
.. _ase-3.6: https://svn.fysik.dtu.dk/projects/ase/tags/3.6.0
.. _ase-3.8.1: https://svn.fysik.dtu.dk/projects/ase/tags/3.8.1


Version 0.8.0
=============

25 May 2011: :trac:`tags/0.8.0 <../tags/0.8.0>`.

* Corresponding ASE release: ase-3.5.1_
* Energy convergence criterion changed from 1 meV/atom to 0.5
  meV/electron.  This was changed in order to allow having no atoms like
  for jellium calculations.
* Linear :ref:`dielectric response <df_theory>` of an extended system
  (RPA and ALDA kernels) can now be calculated.
* :ref:`rpa`.
* Non-selfconsistent calculations with k-points for hybrid functionals.
* Methfessel-Paxton distribution added.
* Text output now shows the distance between planes of grid-points as
  this is what will be close to the grid-spacing parameter *h* also for
  non-orthorhombic cells.
* Exchange-correlation code restructured.  Naming convention for
  explicitely specifying libxc functionals has changed: :ref:`manual_xc`.
* New PAW setups for Rb, Ti, Ba, La, Sr, K, Sc, Ca, Zr and Cs.

.. _ase-3.5.1: https://svn.fysik.dtu.dk/projects/ase/tags/3.5.1


Version 0.7.2
=============

13 August 2010: :trac:`tags/0.7.2 <../tags/0.7.2>`.

* Corresponding ASE release: ase-3.4.1_
* For version 0.7, the default Poisson solver was changed to
  ``PoissonSolver(nn=3)``.  Now, also the Poisson solver's default
  value for ``nn`` has been changed from ``'M'`` to ``3``.

.. _ase-3.4.1:
    https://svn.fysik.dtu.dk/projects/ase/tags/3.4.1

Version 0.7
===========

23 April 2010: :trac:`tags/0.7 <../tags/0.7>`.

* Corresponding ASE release: ase-3.4.0_
* Better and much more efficient handling of non-orthorhombic unit
  cells.  It may actually work now!
* Much better use of ScaLAPACK and BLACS.  All large matrices can now
  be distributed.
* New test coverage pages for all files.
* New default value for Poisson solver stencil: ``PoissonSolver(nn=3)``.
* Much improved MPI module (:ref:`communicators`).
* Self-consistent Meta GGA.
* New :ref:`PAW setup tar-file <setups>` now contains revPBE setups and
  also dzp basis functions.
* New ``$HOME/.gpaw/rc.py`` configuration file.
* License is now GPLv3+.
* New HDF IO-format.
* :ref:`Advanced GPAW Test System <big-test>` Introduced.

.. _ase-3.4.0:
    https://svn.fysik.dtu.dk/projects/ase/tags/3.4.0

Version 0.6
===========

9 October 2009: :trac:`tags/0.6 <../tags/0.6>`.

* Corresponding ASE release: ase-3.2.0_
* Much improved default parameters.
* Using higher order finite-difference stencil for kinetic energy.
* Many many other improvements like: better parallelization, fewer bugs and
  smaller memory footprint.

.. _ase-3.2.0:
    https://svn.fysik.dtu.dk/projects/ase/tags/3.2.0

Version 0.5
===========

1 April 2009: :trac:`tags/0.5 <../tags/0.5>`.

* Corresponding ASE release: ase-3.1.0_
* `new setups added Bi, Br, I, In, Os, Sc, Te; changed Rb setup <https://trac.fysik.dtu.dk/projects/gpaw/changeset/3612>`_.
* `memory estimate feature is back <https://trac.fysik.dtu.dk/projects/gpaw/changeset/3575>`_

.. _ase-3.1.0:
    https://svn.fysik.dtu.dk/projects/ase/tags/3.1.0

Version 0.4
===========

13 November 2008: :trac:`tags/0.4 <../tags/0.4>`.

* Corresponding ASE release: ase-3.0.0_
* Now using ASE-3 and numpy.
* TPSS non self-consistent implementation.
* LCAO mode.
* VdW-functional now coded in C.
* Added atomic orbital basis generation scripts.
* Added an Overlap object, and moved apply_overlap and apply_hamiltonian
  from Kpoint to Overlap and Hamiltonian classes.

* Wannier code much improved.
* Experimental LDA+U code added.
* Now using libxc.
* Many more setups.
* Delta scf calculations.

* Using localized functions will now no longer use MPI group
  communicators and blocking calls to MPI_Reduce and MPI_Bcast.
  Instead non-blocking sends/receives/waits are used.  This will
  reduce syncronization time for large parallel calculations.

* More work on LB94.
* Using LCAO code forinitial guess for grid calculations.
* TDDFT.
* Moved documentation to Sphinx.
* Improved metric for Pulay mixing.
* Porting and optimization for BlueGene/P.
* Experimental Hartwigsen-Goedecker-Hutter pseudopotentials added.
* Transport calculations with LCAO.

.. _ase-3.0.0:
    https://svn.fysik.dtu.dk/projects/ase/tags/3.0.0

Version 0.3
===========

19 December 2007: :trac:`tags/0.3 <../tags/0.3>`.
