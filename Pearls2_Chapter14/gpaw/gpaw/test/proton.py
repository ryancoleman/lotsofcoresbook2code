from __future__ import print_function
from ase import Atoms
from gpaw import GPAW

a = 4.5
H = Atoms('H', [(a / 2, a / 2, a / 2)],
          pbc=0,
          cell=(a, a, a))
calc = GPAW(nbands=1, h=0.2, charge=1)
H.set_calculator(calc)
print(H.get_potential_energy() + calc.get_reference_energy())
assert abs(H.get_potential_energy() + calc.get_reference_energy()) < 0.014
