import numpy as np

from ase.structure import molecule
from ase.io import read, write

a = molecule('CO2')
f = np.array([[1,2,3],[4,5,6],[7,8,9]])
a.set_array('test', f)

write('test.cfg', a)

b = read('test.cfg')
assert np.all(b.get_array('test') == f)
