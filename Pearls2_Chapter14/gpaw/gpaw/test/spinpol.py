from ase import Atom, Atoms
from gpaw import GPAW, FermiDirac
from gpaw.test import equal

a = 4.0
n = 16
hydrogen = Atoms([Atom('H')], cell=(a, a, a), pbc=True)
hydrogen.center()
calc = GPAW(gpts=(n, n, n), nbands=1, convergence={'energy': 1e-5},
            occupations=FermiDirac(0.0))
hydrogen.set_calculator(calc)
e1 = hydrogen.get_potential_energy()
niter1 = calc.get_number_of_iterations()
hydrogen.set_initial_magnetic_moments([1.0])
e2 = hydrogen.get_potential_energy()
niter2 = calc.get_number_of_iterations()
de = e1 - e2
print(de)
equal(de, 0.7871, 1.e-4)

energy_tolerance = 0.00006
niter_tolerance = 0
equal(e1, -0.499854, energy_tolerance)
equal(e2, -1.287, energy_tolerance)
