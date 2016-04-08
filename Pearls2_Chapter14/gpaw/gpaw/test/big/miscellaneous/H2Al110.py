from __future__ import print_function
from ase import Atom, Atoms
from gpaw import GPAW
from gpaw.test import equal

a = 4.00
d = a / 2**0.5
z = 1.1
b = 1.5

slab = Atoms([Atom('Al', (0, 0, 0)),
                    Atom('Al', (a, 0, 0)),
                    Atom('Al', (a/2, d/2, -d/2)),
                    Atom('Al', (3*a/2, d/2, -d/2)),
                    Atom('Al', (0, 0, -d)),
                    Atom('Al', (a, 0, -d)),
                    Atom('Al', (a/2, d/2, -3*d/2)),
                    Atom('Al', (3*a/2, d/2, -3*d/2)),
                    Atom('Al', (0, 0, -2*d)),
                    Atom('Al', (a, 0, -2*d)),
                    Atom('H', (a/2-b/2, 0, z)),
                    Atom('H', (a/2+b/2, 0, z))],
                   cell=(2*a, d, 5*d), pbc=(1, 1, 1))
calc = GPAW(h=0.25, nbands=28, kpts=(2, 6, 1),
            convergence={'eigenstates': 1e-5})
slab.set_calculator(calc)
e = slab.get_potential_energy()
niter = calc.get_number_of_iterations()
assert len(calc.get_k_point_weights()) == 3

for i in range(1):
    slab.positions[-2, 0] -= 0.01
    slab.positions[-1, 0] += 0.01
    e = slab.get_potential_energy()

print(e, niter)
energy_tolerance = 0.00015
niter_tolerance = 0
equal(e, -44.69217, energy_tolerance)
