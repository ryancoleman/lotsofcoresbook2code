from __future__ import print_function
import numpy as np
from ase import Atoms
from gpaw import GPAW
from gpaw.xc.sic import SIC
from gpaw.test import equal
a = 7.0
atom = Atoms('N', magmoms=[3], cell=(a, a, a))
molecule = Atoms('N2', positions=[(0, 0, 0), (0, 0, 1.14)], cell=(a, a, a))
atom.center()
molecule.center()

calc = GPAW(xc=SIC(),
            eigensolver='rmm-diis',
            h=0.17,
            txt='n2.sic.new3b.txt',
            setups='hgh')

atom.set_calculator(calc)
e1 = atom.get_potential_energy()

molecule.set_calculator(calc)
e2 = molecule.get_potential_energy()
F_ac = molecule.get_forces()
print(2 * e1 - e2)
print(F_ac)
