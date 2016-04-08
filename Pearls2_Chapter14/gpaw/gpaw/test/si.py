from __future__ import print_function
from gpaw import GPAW, restart, FermiDirac
from ase import Atoms
from ase.calculators.test import numeric_force
from gpaw.test import equal
import numpy as np

a = 5.404
bulk = Atoms(symbols='Si8',
             scaled_positions=[(0, 0, 0),
                               (0, 0.5, 0.5),
                               (0.5, 0, 0.5),
                               (0.5, 0.5, 0),
                               (0.25, 0.25, 0.25),
                               (0.25, 0.75, 0.75),
                               (0.75, 0.25, 0.75),
                               (0.75, 0.75, 0.25)],
             pbc=True, cell=(a, a, a))
n = 20
calc = GPAW(gpts=(n, n, n),
            nbands=8*3,
            occupations=FermiDirac(width=0.01),
            verbose=1,
            kpts=(1, 1, 1))
bulk.set_calculator(calc)
e1 = bulk.get_potential_energy()
niter1 = calc.get_number_of_iterations()
eigs = calc.get_eigenvalues(kpt=0)
calc.write('temp.gpw')

bulk, calc = restart('temp.gpw', fixdensity=True)
e2 = bulk.get_potential_energy()
eigs2 = calc.get_eigenvalues(kpt=0)
print('Orginal', eigs)
print('Fixdensity', eigs2)
print('Difference', eigs2-eigs)

assert np.fabs(eigs2 - eigs)[:-1].max() < 3e-5
equal(e1, -36.7667, 0.001)
equal(e1, e2, 1e-10)
