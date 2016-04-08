import numpy as np
from ase import Atom, Atoms
from gpaw import GPAW, Mixer, RMM_DIIS
from gpaw.test import equal

a = 4.0
n = 20
d = 1.0
x = d / 3**0.5
atoms = Atoms([Atom('C', (0.0, 0.0, 0.0)),
               Atom('H', (x, x, x)),
               Atom('H', (-x, -x, x)),
               Atom('H', (x, -x, -x)),
               Atom('H', (-x, x, -x))],
              cell=(a, a, a), pbc=True)
calc = GPAW(gpts=(n, n, n), nbands=4, txt='a.txt',
            mixer=Mixer(0.25, 3, 1), eigensolver='rmm-diis')
atoms.set_calculator(calc)
e0 = atoms.get_potential_energy()
niter0 = calc.get_number_of_iterations()

es = RMM_DIIS(blocksize=3)
calc = GPAW(gpts=(n, n, n), nbands=4, txt='b.txt',
            mixer=Mixer(0.25, 3, 1), eigensolver=es)
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()
niter1 = calc.get_number_of_iterations()
equal(e0, e1, 0.000001)
equal(niter0, niter1, 0)
