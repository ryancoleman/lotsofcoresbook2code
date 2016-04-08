from ase.structure import molecule
from gpaw import GPAW

import warnings
# Silence those warnings.
warnings.filterwarnings('ignore', 'Setup for',)

m = molecule('H')
m.center(vacuum=2.0)
calc = GPAW(mode='lcao')
m.set_calculator(calc)
m.get_potential_energy()
calc.write('r.gpw')
calc = GPAW('r.gpw', xc='PBE', idiotproof=False)
calc.get_potential_energy()
