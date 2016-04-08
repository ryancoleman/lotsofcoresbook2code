from __future__ import print_function
from ase import Atoms
from gpaw import GPAW, restart
from gpaw.test import equal

from gpaw.utilities.kspot import AllElectronPotential
if 1:
    be = Atoms(symbols='Be',positions=[(0,0,0)])
    be.center(vacuum=5)
    calc = GPAW(gpts=(64,64,64), xc='LDA', nbands=1) #0.1 required for accuracy
    be.set_calculator(calc)
    e = be.get_potential_energy()
    niter = calc.get_number_of_iterations()
    #calc.write("be.gpw")

    energy_tolerance = 0.00001
    niter_tolerance = 0
    equal(e, 0.00246471, energy_tolerance)

#be, calc = restart("be.gpw")
AllElectronPotential(calc).write_spherical_ks_potentials('bepot.txt')
f = open('bepot.txt')
lines = f.readlines()
f.close()
mmax = 0
for l in lines:
    mmax = max(abs(eval(l.split(' ')[3])), mmax)

print("Max error: ", mmax)
assert mmax<0.009
