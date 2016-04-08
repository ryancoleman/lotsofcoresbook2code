.. _iron_exercise:

====================================
Electron spin and magnetic structure
====================================

As an example of spin polarized calculations, we'll study Fe(bcc) in a
two-atom unit cell, i.e. a simple cubic Bravais lattice with a basis,
where one Fe atom sits at the origin, the other in the middle of the
unit cell. We'll stick to the experimental lattice constant *a* = 2.87
Å.  The atomic term of iron is [Ar]3d\ :sup:`6`\ 4s\ :sup:`2`, i.e. 8
valence electrons/atom is included in the calculation.

* Use Hunds rule (maximum polarization rule) to calculate
  the magnetic moment of an isolated Fe atom.  Draw schematically the
  one-electron eigenvalues for spin up and down on an energy axis,
  along with electron populations.

* We'll make three calculations for bulk iron:

  1) A non-magnetic calculation
  2) A ferro-magnetic calculation (aligned atomic moments)
  3) An anti ferro-magnetic calculation (antiparallel atomic moments).

* How many bands are needed?  Assuming the atoms polarize
  maximally (as the isolated atoms).  For metals, always have at least
  5 extra bands to allow for uneven filling of states for different
  **k**-points.

* One should *help* a magnetic calculation by providing initial
  magnetic moments for the atoms using the ``magmoms`` argument to the
  :class:`~ase.atoms.Atoms` constructor or by using the
  :meth:`~ase.atoms.Atoms.set_initial_magnetic_moments` method of the
  :class:`~ase.atoms.Atoms` object. This
  is necessary to find magnetic states.  Choose the magnetic
  moment close to the expected/desired magnetic state of your system
  (the experimental value is 2.22 per atom). The initial magnetic
  moment is relaxed during the self consistency cycles. 
  Note that for a spin polarized calculation, each iteration step takes 
  twice the time compared to a spin paired calculation.

* For each of the three magnetic phases ferro, antiferro
  and nonmagnetic, write down sensible guesses for initial magnetic
  moment parameters: magnetic moment for each of the two atoms in the
  unit cell.

Start with this script: :download:`ferro.py`.

Compare the energies of the three magnetic phases:

* Experimentally, the ferromagnetic phase is most stable.
  Is this reproduced for LDA and GGA?  Instead of repeating the three
  calculations using PBE, you can estimate the PBE numbers from the LDA
  densities you already have.  This is done in this script:
  :download:`PBE.py`.

* Compare the calculated magnetic moment for the
  ferromagnetic phase with the experimental value.  You can find the
  calculated value in the text output, or by using the
  :meth:`~gpaw.aseinterface.GPAW.get_magnetic_moments` method of the
  calculator object.
