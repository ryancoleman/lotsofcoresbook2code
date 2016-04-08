from ase import Atoms
from gpaw import GPAW

H2 = Atoms('H2', positions=[(0, 0, 0), (1, 0, 0)])
H2.set_cell((3, 3.1, 3.2))
H2.center()
calc = GPAW(convergence={'forces': 1e-8,
                         'density': 100,
                         'energy': 100,
                         'eigenstates': 100})
H2.set_calculator(calc)
H2.get_potential_energy()
assert 17 < calc.iter < 19
