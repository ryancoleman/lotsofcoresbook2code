from __future__ import print_function
from ase import Atom, Atoms
from gpaw import GPAW
from gpaw.test import equal

a = 6.0
calc = GPAW(gpts=(32, 36, 32), nbands=4)
O = Atoms([Atom('O', (a/2, a/2 + 0.5, a/2), magmom=2)],
          pbc=False, cell=(a, a + 1, a), calculator=calc)
e0 = O.get_potential_energy()
niter0 = calc.get_number_of_iterations()

calc.set(charge=1)

e1 = O.get_potential_energy()
niter1 = calc.get_number_of_iterations()

print(e1 - e0)
assert abs(e1 - e0 - 13.989) < 0.04

energy_tolerance = 0.0002
niter_tolerance = 5
equal(e0, -1.88477, energy_tolerance)
equal(e1, 12.11080, energy_tolerance)

# The first ionization energy for LDA oxygen is from this paper:
# In-Ho Lee, Richard M. Martin, Phys. Rev. B 56 7197 (1997)
