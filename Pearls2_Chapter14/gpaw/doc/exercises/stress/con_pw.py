from __future__ import print_function
from ase import Atoms
from ase.lattice import bulk
from gpaw import GPAW
from gpaw import PW

cell = bulk('Si', 'fcc', a=5.421).get_cell()
a = Atoms('Si2', cell=cell, pbc=True,
          scaled_positions=((0, 0, 0), (0.25, 0.25, 0.25)))

for x in [100, 200, 300, 400, 500, 600, 700, 800]:
    # for x in [0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.1]:
    calc = GPAW(mode=PW(x),
                # h=x,
                xc='PBE',
                kpts=(4, 4, 4),
                txt='convergence_%s.txt' % x)

    a.set_calculator(calc)

    print(x, a.get_potential_energy())
