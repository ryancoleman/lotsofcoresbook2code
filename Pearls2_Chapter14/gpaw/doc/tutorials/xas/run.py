from math import pi, cos, sin
from ase import Atoms
from gpaw import GPAW, setup_paths
setup_paths.insert(0, '.')

a = 12.0  # use a large cell

d = 0.9575
t = pi / 180 * 104.51
atoms = Atoms('OH2',
            [(0, 0, 0),
             (d, 0, 0),
             (d * cos(t), d * sin(t), 0)],
            cell=(a, a, a))
atoms.center()
calc = GPAW(nbands=-30,
            h=0.2,
            txt='h2o_xas.txt',
            setups={'O': 'hch1s'})
# the number of unoccupied stated will determine how
# high you will get in energy

atoms.set_calculator(calc)
e = atoms.get_potential_energy()

calc.write('h2o_xas.gpw')
