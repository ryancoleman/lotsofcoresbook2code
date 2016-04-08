from __future__ import print_function
from ase import Atoms
from ase.lattice import bulk
from ase.optimize.bfgs import BFGS
from ase.constraints import UnitCellFilter
from gpaw import GPAW
from gpaw import PW
import numpy as np

cell = bulk('Si', 'fcc', a=6.0).get_cell()
# Experimental Lattice constant is a=5.421 A
a = Atoms('Si2', cell=cell, pbc=True,
          scaled_positions=((0, 0, 0), (0.25, 0.25, 0.25)))

calc = GPAW(xc='PBE',
            mode=PW(400),
            # mode=PW(400, cell=a.get_cell()),  # fix no of planewaves!
            kpts=(4, 4, 4),
            # convergence={'eigenstates': 1.e-10},  # converge tightly!
            txt='stress.txt')
a.set_calculator(calc)

uf = UnitCellFilter(a)
relax = BFGS(uf)
relax.run(fmax=0.05)  # Consider much tighter fmax!

a = np.dot(a.get_cell()[0], a.get_cell()[0])**0.5 * 2**0.5
print('Relaxed lattice parameter: a = %s A' % a)
