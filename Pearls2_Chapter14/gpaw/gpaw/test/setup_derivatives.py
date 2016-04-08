from __future__ import print_function
from ase.structure import molecule
from gpaw import GPAW, PW
from gpaw.atom.derivatives import derivatives

h2 = molecule('H2O')
h2.center(vacuum=3)
h2.calc = GPAW(txt='h2o700.txt', mode=PW(700), setups='new')
print(derivatives(h2))
