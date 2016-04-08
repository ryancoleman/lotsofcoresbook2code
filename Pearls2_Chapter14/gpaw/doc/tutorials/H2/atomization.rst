.. _atomization_energy:

===================================
Calculation of atomization energies
===================================

**Warning**: mainstream DFT is unable to describe correctly the
electronic states of isolated atoms (especially of transition metals).
See http://dx.doi.org/10.1063/1.2723118 .
This usually manifests itself as SCF convergence problems. Please consult
literature before reporting such problems on the mailing lists.

The following script will calculate the atomization
energy of a hydrogen molecule:

.. literalinclude:: atomize.py

First, an ``Atoms`` object containing one hydrogen atom with a
magnetic moment of one, is created.  Next, a GPAW calculator is
created.  The calculator will do a calculation
using the PBE exchange-correlation functional, and output from the
calculation will be written to a file ``H.out``.  The calculator is
hooked on to the ``atom`` object, and the energy is calculated (the
``e1`` variable).  Finally, the result of the calculation
(wavefunctions, densities, ...)  is saved to a file.

The ``molecule`` object is defined, holding the hydrogen molecule at
the experimental lattice constant.  The calculator used for the atom
calculation is used again for the molecule caluclation - only the
filename for the output file needs to be changed to ``H2.out``.  We
extract the energy into the ``e2`` variable.

If we run this script, we get the folowing result:

.. literalinclude:: atomization.txt

According to Blaha *et al.* [1]_, an all-electron calculation with PBE
gives an atomization energy of 4.54 eV, which is in perfect agreement with
our PAW result.

The energy of the spin polarized hydrogen atom is -1.09 eV.  If we do
the calculation for the atom with ``magmom=0`` and ``hund=False``, then we get
almost 0 eV.  This number should converge to exactly zero for a very
large cell and a very high grid-point density, because the energy of a
non spin-polarized hydrogen atom is the reference energy.

.. [1] *S. Kurth, J. P. Perdew, and P. Blaha*, Int. J. Quantum
       Chem. **75**, 889 (1999)
