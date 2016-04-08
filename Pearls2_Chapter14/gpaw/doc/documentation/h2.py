from __future__ import print_function
from ase import Atoms
from gpaw import GPAW

d = 0.74
a = 6.0

atoms = Atoms('H2',
              positions=[(0, 0, 0),
                         (0, 0, d)],
              cell=(a, a, a))
atoms.center()

calc = GPAW(nbands=2, txt='h2.txt')
atoms.set_calculator(calc)
print(atoms.get_forces())
