Getting started with GPAW
=========================

In this exercise we will calculate structures and binding energies for
simple molecules.

Performing a structure optimization
-----------------------------------

A *structure optimization*, also called a *relaxation*, is a series of
calculations used to determine the minimum-energy structure of a given
system.  This involves multiple calculations of the atomic forces
:math:`\mathbf F^a = -\tfrac{\partial E}{\partial \mathbf R^a}` with
respect to the atomic positions :math:`\mathbf R^a` as the atoms
are moved downhill according to an optimization algorithm.

The following script uses the :mod:`EMT calculator <ase.calculators.emt>`
to optimize the structure of :mol:`H_2`.

.. literalinclude:: h2.emt.py

This is the first ASE script we have seen so far, so a few comments
are in order:

* At the top is a series of *import statements*.  These load the
  python modules we are going to use.
* An :class:`~ase.atoms.Atoms` object is created, specifying an initial
  (possibly bad) guess for the atomic positions.
* An :class:`~ase.calculators.emt.EMT` calculator is created.  A
  *calculator* can evaluate
  quantities such as energies and forces on a collection of atoms.
  There are different kinds of calculators, and EMT is a particularly
  simple one.  The calculator is associated with the :class:`~ase.atoms.Atoms`
  object by calling ``atoms.set_calculator(calc)``.
* An :mod:`optimizer <ase.optimize>` is created and
  associated with the
  :class:`~ase.atoms.Atoms` object.  It is also given an optional argument,
  ``trajectory``, which specifies the name of a file into which the
  positions will be saved for each step in the geometry optimization.
* Finally the call ``opt.run(fmax=0.05)`` will run the
  optimization algorithm until all atomic forces are below 0.05 eV per
  Ångström.


**Run the above structure optimization.**

This will print the (decreasing) total energy for each iteration until
it converges, leaving the file :file:`h2.emt.traj` in the working
directory.  Use the command :command:`ase-gui` to view the
trajectory file, showing each step of the optimization.

Structure optimization of :mol:`H_2O` with EMT and GPAW
-------------------------------------------------------

Adapt the above script as needed and calculate the structure of a
:mol:`H_2O` molecule using the EMT calculator.  Note that water is not
a linear molecule.  If you start with a linear molecule, the
minimization may not be able to break the symmetry.  Be sure to
visualize the final configuration to check that it is reasonable.

The empirical EMT potential is fast, but not very accurate for
molecules in particular.  We therefore want to perform this
calculation in GPAW instead.  GPAW uses real-space grids to represent
density and wavefunctions, and the grids exist in a cell.  For this
reason you must set a cell for the :class:`~ase.atoms.Atoms` object.  As a
coarse value let us use a 6 Ångström cell::

  system.set_cell((6.0, 6.0, 6.0))
  system.center()

The cell must be centered in order to prevent atoms from lying too
close to the boundary, as the boundary conditions are zero by default.

Instead of importing and using EMT, we now use GPAW::

  from gpaw import GPAW
  ...
  calc = GPAW()
  ...

Make a copy of your script and adapt it to GPAW, then recalculate the
structure of :mol:`H_2O` (make sure to choose a new filename for the
trajectory file).

During the calculation a lot of text is printed to the terminal.  This
includes the parameters used in the calculation: Atomic positions,
grid spacing, XC functional (GPAW uses LDA by default) and many other
properties.  For each iteration in the self-consistency cycle one line
is printed with the energy and convergence measures.  After the
calculation the energy contributions, band energies and forces are
listed.

Use :command:`ase-gui` to visualize and compare bond lenghts and bond
angles to the EMT result.  Bond lengths and angles are shown
automatically if you select two or three atoms at a time.

Atomization energies
--------------------

Now that we know the structure of :mol:`H_2O`, we can calculate other
interesting properties like the molecule's atomization energy.

The *atomization energy* of a molecule is equal to the total energy of
the molecule minus the sum of the energies of each of its constituent
*isolated* atoms.  For example, the atomization energy of :mol:`H_2` is
:math:`E[\mathrm{H}_2] - 2 E[\mathrm H]`.

GPAW calculations are by default spin-paired, i.e. the spin-up and
spin-down densities are assumed to be equal.  As this is not the case
for isolated atoms, it will be necessary to instruct GPAW to do
something different::

  calc = GPAW(hund=True)

With the ``hund`` keyword, Hund's rule is applied to initialize the
atomic states, and the calculation will be made spin-polarized.

Write a script which calculates the total energy of the isolated O and
H atoms, and calculate the atomization energy of :mol:`H_2O`.

Exchange and correlation functionals
------------------------------------

So far we have been using GPAW's default parameters.  The default
exchange-correlation functional is LDA.  This is not very accurate,
and in particular overestimates atomization energies.  You can specify
different XC functionals to the calculator using
:samp:`GPAW(xc={name})`, where :samp:`{name}` is a string such as
``'LDA'``, ``'PBE'`` or ``'RPBE'``.

Calculate the atomization energy of :mol:`H_2O` with LDA and PBE (just reuse
the geometry from the LDA optimization, i.e. do not repeat the minimization).
