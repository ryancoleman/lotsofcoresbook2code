import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from gpaw import GPAW
from gpaw.wavefunctions.pw import PW
from gpaw.test import equal
from gpaw.mpi import world

a = 2.65
slab = Atoms('Li2',
             [(0, 0, 0), (0, 0, a)],
             cell=(a, a, 3 * a),
             pbc=True)
k = 4
calc = GPAW(mode=PW(200),
            eigensolver='rmm-diis',
            parallel={'band': min(world.size, 4)},
            idiotproof=0,
            kpts=(k, k, 1))
slab.set_calculator(calc)
BFGS(slab).run(fmax=0.01)
assert abs(slab.get_distance(0, 1) - 2.4594) < 0.001, slab.get_distance(0, 1)
