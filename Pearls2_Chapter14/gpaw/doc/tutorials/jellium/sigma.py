from __future__ import print_function
from ase.io import read
e0 = read('bulk.txt').get_potential_energy()
e = read('surface.txt').get_potential_energy()
a = 1.6
L = 10 * a
sigma = (e - L / a * e0) / 2 / a**2
print(('%.2f mev/Ang^2' % (1000 * sigma)))
print(('%.1f erg/cm^2' % (sigma / 6.24150974e-5)))
assert abs(sigma - 0.0054) < 0.0001
