from __future__ import print_function
import os
from ase import Atom, Atoms
from ase.units import Bohr
from gpaw import GPAW
from gpaw.test import equal

a = 7.5 * Bohr
n = 16
atoms = Atoms([Atom('He', (0.0, 0.0, 0.0))], cell=(a, a, a), pbc=True)
calc = GPAW(gpts=(n, n, n), nbands=1, xc='LDA')
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()
niter1 = calc.get_number_of_iterations()
e1ref = calc.get_reference_energy()
de12 = calc.get_xc_difference('PBE')
calc.set(xc='PBE')
e2 = atoms.get_potential_energy()
niter2 = calc.get_number_of_iterations()
e2ref = calc.get_reference_energy()
de21 = calc.get_xc_difference('LDA')
print(e1ref + e1 + de12, e2ref + e2)
print(e1ref + e1, e2ref + e2 + de21)
print(de12, de21)
equal(e1ref + e1 + de12, e2ref + e2, 0.02)
equal(e1ref + e1, e2ref + e2 + de21, 0.025)

calc.write('PBE.gpw')

de21b = GPAW('PBE.gpw').get_xc_difference('LDA')
print(de21, de21b)
equal(de21, de21b, 9e-8)

energy_tolerance = 0.00007
niter_tolerance = 0
equal(e1, -0.0961003634812, energy_tolerance) # svnversion 5252
equal(e2, -0.0790249564625, energy_tolerance) # svnversion 5252
