from __future__ import print_function
"""Test of HGH pseudopotential implementation.

This is the canonical makes-sure-nothing-breaks test, which checks that the
numbers do not change from whatever they were before.

The test runs an HGH calculation on a misconfigured H2O molecule, such that
the forces are nonzero.

Energy is compared to a previous calculation; if it differs significantly,
that is considered an error.

Forces are compared to a previous finite-difference result.
"""

import numpy as np
from ase.structure import molecule
from gpaw import GPAW
from gpaw.utilities import unpack
from gpaw.atom.basis import BasisMaker
from gpaw.test import equal

mol = molecule('H2O')
mol.rattle(0.2)
mol.center(vacuum=2.0)
calc = GPAW(nbands=6,
            gpts=(32, 40, 40),
            setups='hgh',
            convergence=dict(eigenstates=1e-9, density=1e-5, energy=0.3e-5),
            txt='-')
mol.set_calculator(calc)
e = mol.get_potential_energy()
niter = calc.get_number_of_iterations()
F_ac = mol.get_forces()

F_ac_ref = np.array([[ 7.33077718,  3.81069249, -6.07405156],
                     [-0.9079617 , -1.18203514,  3.43777589],
                     [-0.61642527, -0.41889306,  2.332415  ]])

eref = 724.479830538

eerr = abs(e - eref)

print('energy', e)
print('ref energy', eref)
print('energy error', eerr)

print('forces')
print(F_ac)

print('ref forces')
print(F_ac_ref)

ferr = np.abs(F_ac - F_ac_ref).max()
print('max force error', ferr)

fdcheck = False

if fdcheck:
    from ase.calculators.test import numeric_forces
    F_ac_fd = numeric_forces(mol)

    print('Self-consistent forces')
    print(F_ac)
    print('FD forces')
    print(F_ac_fd)
    print()
    print(repr(F_ac_fd))
    print()
    err = np.abs(F_ac - F_ac_fd).max()
    print('max err', err)

wfs = calc.wfs
gd = wfs.gd
psit_nG = wfs.kpt_u[0].psit_nG
dH_asp = calc.hamiltonian.dH_asp

assert eerr < 1e-3, 'energy changed from reference'
assert ferr < 0.015, 'forces do not match FD check'

# Sanity check.  In HGH, the atomic Hamiltonian is constant.
# Also the projectors should be normalized
for a, dH_sp in dH_asp.items():
    dH_p = dH_sp[0]
    K_p = wfs.setups[a].K_p
    #B_ii = wfs.setups[a].B_ii
    #assert np.abs(B_ii.diagonal() - 1).max() < 1e-3
    #print 'B_ii'
    #print wfs.setups[a].B_ii
    # Actually, H2O might not be such a good test, since there'll only
    # be one element in the atomic Hamiltonian for O and zero for H.
    #print 'dH_p', dH_p
    #print 'K_p', K_p

    assert np.abs(dH_p - K_p).max() < 1e-10, 'atomic Hamiltonian changed'

    #h_ii = setup.data.h_ii
    #print 'h_ii', h_ii
    #print 'dH_ii', dH_ii

# Sanity check: HGH is normconserving
for psit_G in psit_nG:
    norm = gd.integrate(psit_G**2) # Around 1e-15 !  Surprisingly good.
    assert abs(1 - norm) < 1e-10, 'Not normconserving'

energy_tolerance = 0.00002
niter_tolerance = 0
equal(e, eref, energy_tolerance) # svnversion 5252
