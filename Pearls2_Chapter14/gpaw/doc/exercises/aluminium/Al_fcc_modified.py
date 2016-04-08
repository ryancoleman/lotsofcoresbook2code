from __future__ import print_function
from numpy.polynomial import Polynomial
import ase.units as u
from ase import Atoms
from gpaw import GPAW, PW

A = [3.9, 4.0, 4.1, 4.2]
E = []
for a in A:
    name = 'bulk-fcc-%.1f' % a
    b = a / 2

    bulk = Atoms('Al',
                 cell=[[0, b, b],
                       [b, 0, b],
                       [b, b, 0]],
                 pbc=True)

    k = 4
    calc = GPAW(mode=PW(300),       # cutoff
                kpts=(k, k, k),     # k-points
                txt=name + '.txt')  # output file

    bulk.set_calculator(calc)

    energy = bulk.get_potential_energy()
    calc.write(name + '.gpw')
    E.append(energy)

p = Polynomial.fit(A, E, 3)
a0 = p.deriv(1).roots()[0]
B = p.deriv(2)(a0) * 4 / 9 / a0 / u.J * u.m**3 * 1e-9  # GPa
print((a0, B))

assert abs(a0 - 3.9924) < 0.001
assert abs(B - 87.22) < 0.1
