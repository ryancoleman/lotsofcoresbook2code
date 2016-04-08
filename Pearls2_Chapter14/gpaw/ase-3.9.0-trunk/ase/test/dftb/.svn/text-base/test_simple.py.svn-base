""" test run for DFTB+ calculator
    the tolerance is extremely loose, beause different sk files
    give different results

"""

from ase.test.dftb import installed
assert installed()

from ase.calculators.dftb import Dftb

import os

from ase.optimize import QuasiNewton

from ase.structure import molecule
test = molecule('H2O')
test.set_calculator(Dftb(label='h2o', atoms=test,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_O='"p"',
                         Hamiltonian_MaxAngularMomentum_H='"s"',
                         ))

dyn = QuasiNewton(test, trajectory='test.traj')
dyn.run(fmax=0.01)
final_energy = test.get_potential_energy()
assert abs(final_energy + 111.141945) < 1.0
files = ['band.out', 'detailed.out', 'dftb_in.hsd', 'dftb_pin.hsd',
         'geo_end.gen', 'geo_end.xyz', 'h2o.out', 'results.tag', 'test.traj',
         'test.traj.bak']
for file in files:
    try:
        os.remove(file)
    except OSError:
        pass
