#!/usr/bin/env python
#

import os

from ase import Atoms
from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.io import write, read

from ase.structure import molecule
test = molecule('H2O')
test.set_calculator(Dftb(label='h2o',atoms=test,
run_manyDftb_steps = True,
Driver_='ConjugateGradient',
Driver_MaxForceComponent='1E-4',
Driver_MaxSteps=1000,
Hamiltonian_MaxAngularMomentum_ = '',
Hamiltonian_MaxAngularMomentum_O = '"p"',
Hamiltonian_MaxAngularMomentum_H = '"s"',
))

dyn = QuasiNewton(test, trajectory='test.traj')
dyn.run(fmax=100, steps=0)
test = read('geo_end.gen')
write('test.final.xyz', test)


