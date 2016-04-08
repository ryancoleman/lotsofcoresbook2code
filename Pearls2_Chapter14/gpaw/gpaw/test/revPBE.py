from ase import Atoms
from gpaw import GPAW
from gpaw.test import equal
from ase.units import Bohr, Hartree

a = 7.5 * Bohr
n = 16
atoms = Atoms('He', [(0.0, 0.0, 0.0)], cell=(a, a, a), pbc=True)
calc = GPAW(gpts=(n, n, n), nbands=1, xc='PBE')
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()
niter1 = calc.get_number_of_iterations()
e1a = calc.get_reference_energy()
calc.set(xc='revPBE')
e2 = atoms.get_potential_energy()
niter2 = calc.get_number_of_iterations()
e2a = calc.get_reference_energy()

equal(e1a, -2.893 * Hartree, 8e-3)
equal(e2a, -2.908 * Hartree, 9e-3)
equal(e1, e2, 4e-3)

energy_tolerance = 0.00005
niter_tolerance = 0
equal(e1, -0.0790449962, energy_tolerance)
equal(e2, -0.08147563, energy_tolerance)
