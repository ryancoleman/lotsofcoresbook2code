import numpy as np

from ase.lattice import bulk
from ase.optimize import BFGS
from ase.io import PickleTrajectory
from ase.constraints import StrainFilter

from gpaw import GPAW, PW


co = bulk('Co')
co.set_initial_magnetic_moments([1.6, 1.6])
co.calc = GPAW(mode=PW(700),
               xc='PBE',
               kpts=(8, 8, 4),
               txt='co.txt')

BFGS(StrainFilter(co)).run(0.005)

a0 = co.cell[0, 0]
c0 = co.cell[2, 2]

traj = PickleTrajectory('co.traj', 'w')
eps = 0.01
for a in a0 * np.linspace(1 - eps, 1 + eps, 3):
    for c in c0 * np.linspace(1 - eps, 1 + eps, 3):
        co.set_cell(bulk('Co', a=a, covera=c / a).cell, scale_atoms=True)
        co.get_potential_energy()
        traj.write(co)
