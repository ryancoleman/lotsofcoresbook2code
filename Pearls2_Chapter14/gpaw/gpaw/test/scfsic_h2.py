import numpy as np
from ase import Atoms
from gpaw import GPAW, restart
from gpaw.xc.sic import SIC
from gpaw.test import equal
a = 6.0
atom = Atoms('H', magmoms=[1.0], cell=(a, a, a))
molecule = Atoms('H2', positions=[(0, 0, 0), (0, 0, 0.737)], cell=(a, a, a))
atom.center()
molecule.center()

calc = GPAW(xc='LDA-PZ-SIC',
            eigensolver='rmm-diis',
            txt='h2.sic.txt',
            setups='hgh')

atom.set_calculator(calc)
e1 = atom.get_potential_energy()

molecule.set_calculator(calc)
e2 = molecule.get_potential_energy()
F_ac = molecule.get_forces()
de = 2 * e1 - e2
#equal(de, 4.5, 0.1)

# Test forces ...

calc.write('H2.gpw', mode='all')
atoms, calc = restart('H2.gpw')
e2b = atoms.get_potential_energy()
equal(e2, e2b, 0.0001)
