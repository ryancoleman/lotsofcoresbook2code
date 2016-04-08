.. module:: ase.constraints
   :synopsis: Constraining some degrees of freedom

===========
Constraints
===========

When performing minimizations or dynamics one may wish to keep some
degrees of freedom in the system fixed. One way of doing this is by
attaching constraint object(s) directly to the atoms object.

Important: setting constraints will freeze the corresponding atom positions.
Changing such atom positions can be achieved:

- by directly setting the :attr:`~ase.atoms.Atoms.positions` attribute
  (see example of setting :ref:`atoms_special_attributes`),

- alternatively, by removing the constraints first::

    del atoms.constraints

  or::

    atoms.set_constraint()

  and using the :meth:`~ase.atoms.Atoms.set_positions` method.


The FixAtoms class
==================

This class is used for fixing some of the atoms.

.. class:: FixAtoms(indices=None, mask=None)

You must supply either the indices of the atoms that should be fixed
or a mask. The mask is a list of booleans, one for each atom, being true
if the atoms should be kept fixed.

For example, to fix the positions of all the Cu atoms in a simulation
with the indices keyword:

>>> c = FixAtoms(indices=[atom.index for atom in atoms if atom.symbol == 'Cu'])
>>> atoms.set_constraint(c)

or with the mask keyword:

>>> c = FixAtoms(mask=[atom.symbol == 'Cu' for atom in atoms])
>>> atoms.set_constraint(c)


The FixBondLength class
=======================

This class is used to fix the distance between two atoms specified by
their indices (*a1* and *a2*)

.. class:: FixBondLength(a1, a2)

Example of use::

  >>> c = FixBondLength(0, 1)
  >>> atoms.set_constraint(c)

In this example the distance between the atoms
with indices 0 and 1 will be fixed in all following dynamics and/or
minimizations performed on the *atoms* object.

This constraint is useful for finding minimum energy barriers for
reactions where the path can be described well by a single bond
length (see the :ref:`mep2` tutorial).

Important: If fixing multiple bond lengths, use the FixBondLengths class
below, particularly if the same atom is fixed to multiple partners.


The FixBondLengths class
========================

More than one bond length can be fixed by using this class. Especially
for cases in which more than one bond length constraint is applied on
the same atom. It is done by specifying the indices of the two atoms
forming the bond in pairs.

.. class:: FixBondLengths(pairs)

Example of use::

  >>> c = FixBondLengths([[0, 1], [0, 2]])
  >>> atoms.set_constraint(c)

Here the distances between atoms with indices 0 and 1 and atoms with
indices 0 and 2 will be fixed. The constraint is for the same purpose
as the FixBondLength class.


The FixedLine class
====================

.. autoclass:: FixedLine


The FixedPlane class
====================

.. autoclass:: FixedPlane

Example of use: :ref:`constraints_diffusion_tutorial`.


The FixedMode class
===================

.. autoclass:: FixedMode

A mode is a list of vectors specifying a direction for each atom. It often
comes from :meth:`ase.vibrations.Vibrations.get_mode`.


The Hookean class
=================

This class of constraints, based on Hooke's Law, is generally used to
conserve molecular identity in optimization schemes and can be used in three
different ways. In the first, it applies a Hookean restorative force between
two atoms if the distance between them exceeds a threshold. This is useful to
maintain the identity of molecules in quenched molecular dynamics, without
changing the degrees of freedom or violating conservation of energy. When the
distance between the two atoms is less than the threshold length, this
constraint is completely inactive.

The below example tethers atoms at indices 3 and 4 together::

  >>> c = Hookean(a1=3, a2=4, rt=1.79, k=5.)
  >>> atoms.set_constraint(c)

Alternatively, this constraint can tether a single atom to a point in space,
for example to prevent the top layer of a slab from subliming during a
high-temperature MD simulation. An example of tethering atom at index 3 to its
original position:

>>> c = Hookean(a1=3, a2=atoms[3].position, rt=0.94, k=2.)
>>> atoms.set_constraint(c)

Reasonable values of the threshold (rt) and spring constant (k) for some
common bonds are below.

.. list-table::

  * - Bond
    - rt (Angstroms)
    - k (eV Angstrom^-2)
  * - O-H
    - 1.40
    - 5
  * - C-O
    - 1.79
    - 5
  * - C-H
    - 1.59
    - 7
  * - C=O
    - 1.58
    - 10
  * - Pt sublimation
    - 0.94
    - 2
  * - Cu sublimation
    - 0.97
    - 2

A third way this constraint can be applied is to apply a restorative force if an atom crosses a plane in space. For example::

  >>> c = Hookean(a1=3, a2=(0, 0, 1, -7), k=10.)
  >>> atoms.set_constraint(c)

This will apply a restorative force on atom 3 in the downward direction of magnitude k * (atom.z - 7) if the atom's vertical position exceeds 7 Angstroms. In other words, if the atom crosses to the (positive) normal side of the plane, the force is applied and directed towards the plane. (The same plane with the normal direction pointing in the -z direction would be given by (0, 0, -1, 7).)

For an example of use, see the :ref:`mhtutorial` tutorial.

.. note::

  In previous versions of ASE, this was known as the BondSpring constraint.

  
The FixInternals class
======================

This class allows to fix an arbitrary number of bond lengths, angles
and dihedral angles. The defined constraints are satisfied self
consistently. To define the constraints one needs to specify the
atoms object on which the constraint works (needed for atomic
masses), a list of bond, angle and dihedral constraints.
Those constraint definitions are always list objects containing
the value to be set and a list of atomic indices. The epsilon value
specifies the accuracy to which the constraints are fulfilled.

.. class:: FixInternals(atoms, bonds=[bond1, bond2], \
    angles=[angle1], dihedrals=[dihedral1, dihedral2], epsilon=1.e-7)

Example of use::

  >>> bond1 = [1.20, [1, 2]]
  >>> angle_indices1 = [2, 3, 4]
  >>> dihedral_indices1 = [2, 3, 4, 5]
  >>> angle1 = [atoms.get_angle(angle_indices1), angle_indices1]
  >>> dihedral1 = [atoms.get_dihedral(dihedral_indices1), \
    dihedral_indices1]
  >>> c = FixInternals(atoms, bonds=[bonds1], angles=[angles1], \
    dihedrals=[dihedral1])
  >>> atoms.set_constraint(c)

This example defines a bond an angle and a dihedral angle constraint
to be fixed at the same time.


Combining constraints
=====================

It is possible to supply several constraints on an atoms object. For
example one may wish to keep the distance between two nitrogen atoms
fixed while relaxing it on a fixed ruthenium surface::

  >>> pos = [[0.00000, 0.00000,  9.17625],
  ...        [0.00000, 0.00000, 10.27625],
  ...        [1.37715, 0.79510,  5.00000],
  ...        [0.00000, 3.18039,  5.00000],
  ...        [0.00000, 0.00000,  7.17625],
  ...        [1.37715, 2.38529,  7.17625]]
  >>> unitcell = [5.5086, 4.7706, 15.27625]

  >>> atoms = Atoms(positions=pos,
  ...               symbols='N2Ru4',
  ...               cell=unitcell,
  ...               pbc=[True,True,False])

  >>> fa = FixAtoms(mask=[a.symbol == 'Ru' for a in atoms])
  >>> fb = FixBondLength(0, 1)
  >>> atoms.set_constraint([fa, fb])

When applying more than one constraint they are passed as a list in
the :meth:`~ase.atoms.Atoms.set_constraint` method, and they will be applied
one after the other.

Important: If wanting to fix the length of more than one bond in the
simulation, do not supply a list of :class:`FixBondLength`
instances; instead, use a single instance of
:class:`FixBondLengths`.


Making your own constraint class
================================

A constraint class must have these two methods:

.. method:: adjust_positions(oldpositions, newpositions)

   Adjust the *newpositions* array inplace.

.. method:: adjust_forces(positions, forces)

   Adjust the *forces* array inplace.


A simple example::

  import numpy as np
  class MyConstraint:
      """Constrain an atom to move along a given direction only."""
      def __init__(self, a, direction):
          self.a = a
          self.dir = direction / sqrt(np.dot(direction, direction))

      def adjust_positions(self, oldpositions, newpositions):
          step = newpositions[self.a] - oldpositions[self.a]
          step = np.dot(step, self.dir)
          newpositions[self.a] = oldpositions[self.a] + step * self.dir

      def adjust_forces(self, positions, forces):
          forces[self.a] = self.dir * np.dot(forces[self.a], self.dir)

A constraint can optionally have two additional methods, which
will be ignored if missing:

.. method:: adjust_momenta(positions, momenta)

   Adjust the *momenta* array inplace.

.. method:: adjust_potential_energy(positions, energy)

   Provide the difference in the *potential energy* due to the constraint.
   (Note that inplace adjustment is not possible for energy, which is a
   float.)


The Filter class
================

Constraints can also be applied via filters, which acts as a wrapper
around an atoms object. A typical use case will look like this::

   -------       --------       ----------
  |       |     |        |     |          |
  | Atoms |<----| Filter |<----| Dynamics |
  |       |     |        |     |          |
   -------       --------       ----------

and in Python this would be::

  >>> atoms = Atoms(...)
  >>> filter = Filter(atoms, ...)
  >>> dyn = Dynamics(filter, ...)


This class hides some of the atoms in an Atoms object.

.. class:: Filter(atoms, indices=None, mask=None)

You must supply either the indices of the atoms that should be kept
visible or a mask. The mask is a list of booleans, one for each atom,
being true if the atom should be kept visible.

Example of use::

  >>> from ase import Atoms, Filter
  >>> atoms=Atoms(positions=[[ 0    , 0    , 0],
  ...                        [ 0.773, 0.600, 0],
  ...                        [-0.773, 0.600, 0]],
  ...             symbols='OH2')
  >>> f1 = Filter(atoms, indices=[1, 2])
  >>> f2 = Filter(atoms, mask=[0, 1, 1])
  >>> f3 = Filter(atoms, mask=[a.Z == 1 for a in atoms])
  >>> f1.get_positions()
  [[ 0.773  0.6    0.   ]
   [-0.773  0.6    0.   ]]

In all three filters only the hydrogen atoms are made
visible.  When asking for the positions only the positions of the
hydrogen atoms are returned.


The UnitCellFilter class
========================

.. autoclass:: UnitCellFilter

The StrainFilter class
======================

.. autoclass:: StrainFilter
