from __future__ import print_function
from sys import argv
import matplotlib.pyplot as plt

from ase.dft import STM
from gpaw import restart

filename = argv[1]
z0 = 8
bias = 1.0

atoms, calc = restart(filename, txt=None)

stm = STM(atoms, symmetries=[0, 1, 2])
c = stm.get_averaged_current(bias, z0)

print('Average current at z=%f: %f' % (z0, c))

# Get 2d array of constant current heights:
x, y, h = stm.scan(bias, c)

print('Min: %.2f Ang, Max: %.2f Ang' % (h.min(), h.max()))

plt.contourf(x, y, h, 40)
plt.hot()
plt.colorbar()
plt.show()
