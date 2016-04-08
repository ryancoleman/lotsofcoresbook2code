"""Make sure we get an exception when an atom is too close to the boundary."""
from ase import Atoms
from gpaw import GPAW
from gpaw.grid_descriptor import GridBoundsError
a = 4.0
x = 0.1
hydrogen = Atoms('H', [(x, x, x)],
                 cell=(a, a, a),
                 calculator=GPAW(maxiter=7))
try:
    e1 = hydrogen.get_potential_energy()
except GridBoundsError:
    pass
else:
    assert False
