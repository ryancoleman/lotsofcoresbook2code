from __future__ import print_function
from ase import Atom, Atoms
from ase.io import read
from gpaw import GPAW, FermiDirac
from gpaw.test import equal
a = 2.0
calc = GPAW(gpts=(12, 12, 12), txt='H.txt', occupations=FermiDirac(0.0))
H = Atoms('H',
          cell=(a, a, a),
          pbc=True,
          calculator=calc)
e0 = H.get_potential_energy()

H = read('H.txt')
equal(H.get_potential_energy(), e0, 1e-6)

energy_tolerance = 0.00007
equal(e0, -6.55685, energy_tolerance)

print(calc.get_xc_functional())
