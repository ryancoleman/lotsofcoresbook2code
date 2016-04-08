import numpy as np
from ase.structure import molecule
from ase.dft import Wannier
from gpaw import GPAW
from gpaw.test import equal

# Test of ase wannier using gpaw

calc = GPAW(gpts=(32, 32, 32), nbands=4)
atoms = molecule('H2', calculator=calc)
atoms.center(vacuum=3.)
e = atoms.get_potential_energy()
niter = calc.get_number_of_iterations()

pos = atoms.positions + np.array([[0, 0, .2339], [0, 0, -.2339]])
com = atoms.get_center_of_mass()

wan = Wannier(nwannier=2, calc=calc, initialwannier='bloch')
equal(wan.get_functional_value(), 2.964, 1e-3)
equal(np.linalg.norm(wan.get_centers() - [com, com]), 0, 1e-4)

wan = Wannier(nwannier=2, calc=calc, initialwannier='projectors')
equal(wan.get_functional_value(), 3.100, 2e-3)
equal(np.linalg.norm(wan.get_centers() - pos), 0, 1e-3)

wan = Wannier(nwannier=2, calc=calc, initialwannier=[[0, 0, .5], [1, 0, .5]])
equal(wan.get_functional_value(), 3.100, 2e-3)
equal(np.linalg.norm(wan.get_centers() - pos), 0, 1e-3)

wan.localize()
equal(wan.get_functional_value(), 3.100, 2e-3)
equal(np.linalg.norm(wan.get_centers() - pos), 0, 1e-3)
equal(np.linalg.norm(wan.get_radii() - 1.2393), 0, 2e-3)
eig = np.sort(np.linalg.eigvals(wan.get_hamiltonian().real))
equal(np.linalg.norm(eig - calc.get_eigenvalues()[:2]), 0, 1e-4)

energy_tolerance = 0.00005
niter_tolerance = 0
equal(e, -6.65064, energy_tolerance)
