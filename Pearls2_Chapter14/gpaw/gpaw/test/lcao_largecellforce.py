from __future__ import print_function
# This test calculates the force on the atoms in a hydrogen chain, comparing
# the results to finite-difference values.
#
# The reference result is obtained from an FD calculation, which can be rerun
# by setting the fd boolean below.
#
# The purpose is to test domain decomposition with large cells.  The basis
# functions of one atom are defined to not overlap the rightmost domain
# for z-decompositions of two or more slices.  This will change the set
# of atomic indices stored in BasisFunctions objects and other things
#
# This test also ensures that lcao forces are tested with non-pbc.
import numpy as np
from numpy import array
from ase import Atoms
from ase.units import Bohr
from gpaw import GPAW
from gpaw.atom.basis import BasisMaker

hbasis = BasisMaker('H').generate(1, 0, energysplit=1.8, tailnorm=0.03**.5)
basis = {'H' : hbasis}

atom = Atoms('H')
atom.center(vacuum=.8)
system = atom.repeat((1, 1, 4))

rc = hbasis.bf_j[0].rc

system.center(vacuum=2.5)

calc = GPAW(h=0.23,
            mode='lcao',
            basis=basis,
            convergence={'density':1e-4, 'energy': 1e-7},
            )

system.set_calculator(calc)

F_ac = system.get_forces()

# Check that rightmost domain is in fact outside range of basis functions
from gpaw.mpi import rank, size
if rank == 0 and size > 1:
    assert len(calc.wfs.basis_functions.atom_indices) < len(system)

fd = 0

# Values taken from FD calculation below
# (Symmetry means only z-component may be nonzero)
ref = array([[0.0, 0.0,  4.61734874],
             [0.0, 0.0, -2.74398046],
             [0.0, 0.0,  2.74398027],
             [0.0, 0.0, -4.61734856]])

if fd:
    from ase.calculators.test import numeric_forces
    ref = numeric_forces(system, axes=[2], d=0.002)
    print('Calced')
    print(F_ac)
    print('FD')
    print(ref)
    print(repr(ref))

err = np.abs(F_ac - ref).max()
print('Ref')
print(ref)
print('Calculated')
print(F_ac)
print('Max error', err)
assert err < 6e-4
