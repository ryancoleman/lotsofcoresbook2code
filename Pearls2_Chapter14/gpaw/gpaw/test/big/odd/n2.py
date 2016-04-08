from ase.structure import molecule
from gpaw import GPAW
from gpaw.odd.sic import SIC
from gpaw import extra_parameters
from gpaw.test import equal

extra_parameters['sic'] = True

n = molecule('N')
n.center(vacuum=3)
calc = GPAW(xc=SIC(nspins=2),
            h=0.2,
            convergence=dict(eigenvalues=1e-5, density=1e-3),
            txt='N.txt')
n.set_calculator(calc)
e1 = n.get_potential_energy()

n2 = molecule('N2')
n2.center(vacuum=3)
calc = GPAW(xc=SIC(nspins=1),
            h=0.2,
            convergence=dict(eigenvalues=1e-5, density=1e-3),
            txt='N2.txt')
n2.set_calculator(calc)
e2 = n2.get_potential_energy()

equal(e2 - 2 * e1, -4.5, 0.1)
