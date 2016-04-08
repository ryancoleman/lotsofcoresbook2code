from ase import Atom, Atoms
from gpaw import GPAW
from gpaw.test import equal

a = 4.05
d = a / 2**0.5
bulk = Atoms([Atom('Al', (0, 0, 0)),
              Atom('Al', (0.5, 0.5, 0.5))], pbc=True)
bulk.set_cell((d, d, a), scale_atoms=True)
h = 0.25
calc = GPAW(mode='pw',
            nbands=2*8,
            kpts=(2, 2, 2),
            convergence={'eigenstates': 7.2e-9, 'energy': 1e-5})
bulk.set_calculator(calc)
e0 = bulk.get_potential_energy()
niter0 = calc.get_number_of_iterations()
calc = GPAW(mode='pw',
            nbands=2*8,
            kpts=(2, 2, 2),
            convergence={'eigenstates': 7.2e-9,
                         'energy': 1e-5,
                         'bands': 5 },
            eigensolver='dav')
bulk.set_calculator(calc)
e1 = bulk.get_potential_energy()
niter1 = calc.get_number_of_iterations()
equal(e0, e1, 5.0e-6)

energy_tolerance = 0.00004
niter_tolerance = 0
equal(e0, -6.97798, energy_tolerance)
assert 10 <= niter0 <= 14, niter0
equal(e1, -6.97798, energy_tolerance)
assert 10 <= niter1 <= 24, niter1
