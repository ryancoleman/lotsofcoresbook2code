from __future__ import print_function
import numpy as np

from gpaw.grid_descriptor import GridDescriptor
from gpaw.spline import Spline
import gpaw.mpi as mpi
from gpaw.wavefunctions.pw import PWDescriptor, PWLFC


x = 2.0
rc = 3.5
r = np.linspace(0, rc, 100)

n = 40
a = 8.0
cell_cv = np.array([[a, 0.5, -1], [0, a, 2], [-1, 0, a + 1]])
gd = GridDescriptor((n, n, n), cell_cv, comm=mpi.serial_comm)

a_R = gd.empty()
z = np.linspace(0, n, n, endpoint=False)
a_R[:] = 2 + np.sin(2 * np.pi * z / n)

spos_ac = np.array([(0.15, 0.45, 0.95)])

pd = PWDescriptor(45, gd)
a_G = pd.fft(a_R)

s = Spline(0, rc, 2 * x**1.5 / np.pi * np.exp(-x * r**2))
p = Spline(1, rc, 2 * x**1.5 / np.pi * np.exp(-x * r**2))
d = Spline(2, rc, 2 * x**1.5 / np.pi * np.exp(-x * r**2))

lfc = PWLFC([[s, p, d]], pd)
lfc.set_positions(spos_ac)
b_LG = pd.zeros(9)
lfc.add(b_LG, {0: np.eye(9)})
e1 = pd.integrate(a_G, b_LG)
assert abs(lfc.integrate(a_G)[0] - e1).max() < 1e-11

s1 = []
for i in range(9):
    x = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    x[i] = 1
    s1.append(lfc.stress_tensor_contribution(a_G, {0: x}) -
              np.eye(3) * e1[i])

x = 1e-6
for dist in [[[x, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[0, 0, 0], [0, x, 0], [0, 0, 0]],
             [[0, 0, 0], [0, 0, 0], [0, 0, x]],
             [[0, x, 0], [x, 0, 0], [0, 0, 0]],
             [[0, 0, x], [0, 0, 0], [x, 0, 0]],
             [[0, 0, 0], [0, 0, x], [0, x, 0]]]:
    e = dist + np.eye(3)
    c_cv = np.dot(cell_cv, e)
    gd = GridDescriptor((n, n, n), c_cv, comm=mpi.serial_comm)
    pd = PWDescriptor(45, gd)
    aa_G = a_G / np.linalg.det(e)
    lfc = PWLFC([[s, p, d]], pd)
    lfc.set_positions(spos_ac)
    b_LG = pd.zeros(9)
    lfc.add(b_LG, {0: np.eye(9)})
    e2 = pd.integrate(aa_G, b_LG)
    s2 = (e2 - e1) / x
    error = (np.array(s1) * dist).sum(1).sum(1) / x - s2
    print(s2, abs(error).max())
    assert abs(error).max() < 2e-6
