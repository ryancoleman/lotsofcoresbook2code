from ase.lattice import bulk
from gpaw import GPAW
si = bulk('Si', 'diamond', a=5.5, cubic=not True)
si.set_calculator(GPAW(setups='ah', kpts=(2, 2, 2)))
si.get_potential_energy()
si.calc.write('Si.gpw', 'all')
GPAW('Si.gpw')
