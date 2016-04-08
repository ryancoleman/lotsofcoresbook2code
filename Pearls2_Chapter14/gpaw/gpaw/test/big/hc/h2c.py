"""Hydrogen chain test.

See:

    S. Suhai: *Electron correlation in extended systems: Fourth-order
    many-body perturbation theory and density-functional methods
    applied to an infinite chain of hydrogen atoms*, Phys. Rev. B 50,
    14791-14801 (1994)
    
"""

import numpy as np
from ase import Atoms

from gpaw import GPAW, FermiDirac

k = 8
a = 6.0
d0 = 1.0

h2 = Atoms('H2', cell=(2 * d0, a, a), pbc=True,
           positions=[(0, 0, 0), (d0 + 0.1, 0, 0)])
h2.center(axis=1)
h2.center(axis=2)

h2.calc = GPAW(kpts=(k, 1, 1),
               occupations=FermiDirac(0.01),
               txt='h2c.txt')

for d in np.linspace(0.0, 0.4, 21):
    h2[1].x = d0 + d
    e = h2.get_potential_energy()
    h2.calc.write('h2c-%.2f.gpw' % d)
