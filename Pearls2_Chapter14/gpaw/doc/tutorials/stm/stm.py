from ase.dft.stm import STM
from gpaw import GPAW
calc = GPAW('al111.gpw')
atoms = calc.get_atoms()
stm = STM(atoms)
z = 8.0
bias = 1.0
c = stm.get_averaged_current(bias, z)
x, y, h = stm.scan(bias, c, repeat=(3, 5))
import matplotlib.pyplot as plt
import numpy as np
plt.gca(aspect='equal')
plt.contourf(x, y, h, 40)
plt.hot()
plt.colorbar()
plt.savefig('2d.png')
plt.figure()
a = atoms.cell[0, 0]
x, y = stm.linescan(bias, c, [0, 0], [2 * a, 0])
plt.plot(x, y)
plt.savefig('line.png')
