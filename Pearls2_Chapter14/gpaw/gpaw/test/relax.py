from __future__ import print_function
import os
from ase import Atom, Atoms
from ase.optimize import BFGS
from ase.io import read
from gpaw import GPAW
from gpaw.test import equal

a = 4.    # Size of unit cell (Angstrom)
c = a / 2
d = 0.74  # Experimental bond length
molecule = Atoms('H2',
                 [(c - d / 2, c, c),
                  (c + d / 2, c, c)],
                 cell=(a, a, a),
                 pbc=False)
calc = GPAW(h=0.2, nbands=1, xc='PBE', txt=None)
molecule.set_calculator(calc)
e1 = molecule.get_potential_energy()
niter1 = calc.get_number_of_iterations()
calc.write('H2.gpw')
calc.write('H2a.gpw', mode='all')
molecule.get_forces()
calc.write('H2f.gpw')
calc.write('H2fa.gpw', mode='all')

from time import time
def timer(func, *args, **kwargs):
    t0 = time()
    ret = func(*args, **kwargs)
    return ret, time()-t0

molecule = GPAW('H2.gpw', txt=None).get_atoms()
f1, t1 = timer(molecule.get_forces)
molecule = GPAW('H2a.gpw', txt=None).get_atoms()
f2, t2 = timer(molecule.get_forces)
molecule = GPAW('H2f.gpw', txt=None).get_atoms()
f3, t3 = timer(molecule.get_forces)
molecule = GPAW('H2fa.gpw', txt=None).get_atoms()
f4, t4 = timer(molecule.get_forces)
print('timing:', t1, t2, t3, t4)
assert t2 < 0.6 * t1
assert t3 < 0.5
assert t4 < 0.5
print(f1)
print(f2)
print(f3)
print(f4)
assert sum((f1 - f4).ravel()**2) < 1e-6
assert sum((f2 - f4).ravel()**2) < 1e-6
assert sum((f3 - f4).ravel()**2) < 1e-6

positions = molecule.get_positions()
#                 x-coordinate      x-coordinate
#                 v                 v
d0 = positions[1, 0] - positions[0, 0]
#              ^                 ^
#              second atom       first atom

print('experimental bond length:')
print('hydrogen molecule energy: %7.3f eV' % e1)
print('bondlength              : %7.3f Ang' % d0)

# Find the theoretical bond length:
relax = BFGS(molecule)
relax.run(fmax=0.05)

e2 = molecule.get_potential_energy()
niter2 = calc.get_number_of_iterations()

positions = molecule.get_positions()
#                 x-coordinate      x-coordinate
#                 v                 v
d0 = positions[1, 0] - positions[0, 0]
#              ^                 ^
#              second atom       first atom

print('PBE energy minimum:')
print('hydrogen molecule energy: %7.3f eV' % e2)
print('bondlength              : %7.3f Ang' % d0)


molecule = GPAW('H2fa.gpw', txt='H2.txt').get_atoms()
relax = BFGS(molecule)
relax.run(fmax=0.05)
e2q = molecule.get_potential_energy()
niter2q = calc.get_number_of_iterations()
positions = molecule.get_positions()
d0q = positions[1, 0] - positions[0, 0]
assert abs(e2 - e2q) < 2e-6
assert abs(d0q - d0) < 4e-4

f0 = molecule.get_forces()
del relax, molecule

from gpaw.mpi import world
world.barrier()  # syncronize before reading text output file
f = read('H2.txt').get_forces()
assert abs(f - f0).max() < 5e-6  # 5 digits in txt file

energy_tolerance = 0.00005
niter_tolerance = 0
equal(e1, -6.287873, energy_tolerance)
equal(e2, -6.290744, energy_tolerance)
equal(e2q, -6.290744, energy_tolerance)
