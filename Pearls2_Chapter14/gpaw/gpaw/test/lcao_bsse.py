from __future__ import print_function
from ase.structure import molecule
from gpaw import GPAW
from gpaw.poisson import PoissonSolver
from gpaw.atom.basis import BasisMaker
from gpaw.test import equal

# Tests basis set super position error correction

# Compares a single hydrogen atom to a system of one hydrogen atom
# and one ghost hydrogen atom.  The systems should have identical properties,
# i.e. the ghost orbital should have a coefficient of 0.

b = BasisMaker('H').generate(1, 0, energysplit=0.005)

system = molecule('H2')
system.center(vacuum=6.0)

def prepare(setups):
    calc = GPAW(basis={'H' : b}, mode='lcao',
                setups=setups, h=0.2,
                poissonsolver=PoissonSolver('M', relax='GS', eps=1e-5),
                spinpol=False,
                nbands=1)
    system.set_calculator(calc)
    return calc

calc = prepare({0 : 'paw', 1 : 'ghost'})
system.set_calculator(calc)
e_bsse = system.get_potential_energy()
niter_bsse = calc.get_number_of_iterations()

c_nM = calc.wfs.kpt_u[0].C_nM
print('coefs')
print(c_nM)
print('energy', e_bsse)

# Reference system which is just a hydrogen
sys0 = system[0:1].copy()
calc = prepare('paw')
sys0.set_calculator(calc)
e0 = sys0.get_potential_energy()
niter0 = calc.get_number_of_iterations()
print('e0, e_bsse = ', e0, e_bsse)

# One coefficient should be very small (0.012), the other very large (0.99)
assert abs(1.0 - abs(c_nM[0, 0])) < 0.02
assert abs(c_nM[0, 1]) < 0.02
assert abs(e_bsse - e0) < 2e-3

energy_tolerance = 0.0002
niter_tolerance = 1
equal(e_bsse, 0.0287208853911, energy_tolerance) # svnversion 5252
equal(niter_bsse, 7, niter_tolerance) # svnversion 5252
equal(e0, 0.0299220702846, energy_tolerance) # svnversion 5252
equal(niter0, 6, niter_tolerance) # svnversion 5252
