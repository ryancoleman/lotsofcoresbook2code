from __future__ import print_function
from ase import Atom, Atoms
from gpaw import GPAW, PoissonSolver
from gpaw.test import equal
from gpaw.xc.hybrid import HybridXC

a = 6.  # Size of unit cell (Angstrom)
c = a / 2
# Hydrogen atom:
atom = Atoms([Atom('H', (c, c, c), magmom=1)],
                   cell=(a, a, a), pbc=False)

# gpaw calculator:
calc = GPAW(gpts=(32, 32, 32), nbands=1, xc='PBE', txt='H.txt',
            convergence=dict(eigenstates=3.3e-8))
atom.set_calculator(calc)

e1 = atom.get_potential_energy()
niter1 = calc.get_number_of_iterations()
de1t = calc.get_xc_difference('TPSS')
de1m = calc.get_xc_difference('M06L')
de1x = calc.get_xc_difference(HybridXC('EXX', finegrid=True))
de1xb = calc.get_xc_difference(HybridXC('EXX', finegrid=False))

# Hydrogen molecule:
d = 0.74  # Experimental bond length
molecule = Atoms([Atom('H', (c - d / 2, c, c)),
                  Atom('H', (c + d / 2, c, c))],
                 cell=(a, a, a), pbc=False)

calc.set(txt='H2.txt')
molecule.set_calculator(calc)
e2 = molecule.get_potential_energy()
niter2 = calc.get_number_of_iterations()
de2t = calc.get_xc_difference('TPSS')
de2m = calc.get_xc_difference('M06L')
de2x = calc.get_xc_difference(HybridXC('EXX', finegrid=True))
de2xb = calc.get_xc_difference(HybridXC('EXX', finegrid=False))

print('hydrogen atom energy:     %5.2f eV' % e1)
print('hydrogen molecule energy: %5.2f eV' % e2)
print('atomization energy:       %5.2f eV' % (2 * e1 - e2))
print('atomization energy  TPSS: %5.2f eV' % (2 * (e1 + de1t) - (e2 + de2t)))
print('atomization energy  M06L: %5.2f eV' % (2 * (e1 + de1m) - (e2 + de2m)))
print('atomization energy   EXX: %5.2f eV' % (2 * (e1 + de1x) - (e2 + de2x)))
print('atomization energy   EXX: %5.2f eV' % (2 * (e1 + de1xb) - (e2 + de2xb)))
PBETPSSdifference = (2 * e1 - e2) - (2 * (e1 + de1t) - (e2 + de2t))
PBEM06Ldifference = (2 * e1 - e2) - (2 * (e1 + de1m) - (e2 + de2m))
PBEEXXdifference = (2 * e1 - e2) - (2 * (e1 + de1x) - (e2 + de2x))
PBEEXXbdifference = (2 * e1 - e2) - (2 * (e1 + de1xb) - (e2 + de2xb))
print(PBETPSSdifference)
print(PBEM06Ldifference)
print(PBEEXXdifference)
print(PBEEXXbdifference)
# TPSS value is from JCP 120 (15) 6898, 2004
# e.g. Table VII: DE(PBE - TPSS) = (104.6-112.9)*kcal/mol
# EXX value is from PRL 77, 3865 (1996)
equal(PBETPSSdifference, -0.3599, 0.04)
equal(PBEM06Ldifference, -0.169, 0.01)
equal(PBEEXXdifference, 0.91, 0.005)
equal(PBEEXXbdifference, 0.91, 0.005)

energy_tolerance = 0.0002
niter_tolerance = 0
equal(e1, -1.081638, energy_tolerance)
equal(e2, -6.726356, energy_tolerance)
