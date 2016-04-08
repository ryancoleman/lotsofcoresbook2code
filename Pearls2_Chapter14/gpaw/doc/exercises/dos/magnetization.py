from ase.io import write
from gpaw import GPAW
calc = GPAW('anti.gpw')
atoms = calc.get_atoms()
up = calc.get_pseudo_density(0)
down = calc.get_pseudo_density(1)
zeta = (up - down) / (up + down)
write('magnetization.cube', atoms, data=zeta)
