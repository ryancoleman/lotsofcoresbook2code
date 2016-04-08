from ase import Atoms
from gpaw import GPAW

n = Atoms('N', magmoms=[3])
n.center(vacuum=3.5)

# Calculation with no +U correction:
n.calc = GPAW(mode='lcao',
              basis='dzp',
              txt='no_u.txt',
              xc='PBE')
e1 = n.get_potential_energy()

# Calculation with a correction U=6 eV normalized:
n.calc.set(setups={'N': ':p,6.0'}, txt='normalized_u.txt')
e2 = n.get_potential_energy()

# Calculation with a correction U=6 eV not normalized:
n.calc.set(setups={'N': ':p,6.0,0'}, txt='not_normalized_u.txt')
e3 = n.get_potential_energy()
