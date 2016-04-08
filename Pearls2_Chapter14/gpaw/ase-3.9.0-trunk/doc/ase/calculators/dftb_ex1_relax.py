from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.io import write

from ase.structure import molecule
test = molecule('H2O')
test.set_calculator(Dftb(label='h2o',
                         atoms=test,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_O='"p"',
                         Hamiltonian_MaxAngularMomentum_H='"s"',
                         ))

dyn = QuasiNewton(test, trajectory='test.traj')
dyn.run(fmax=0.01)
write('test.final.xyz', test)



