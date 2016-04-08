from __future__ import print_function
import numpy as np
from gpaw.lfc import LocalizedFunctionsCollection as LFC
from gpaw.grid_descriptor import GridDescriptor
from gpaw.spline import Spline
gd = GridDescriptor([20, 16, 16], [(4, 2, 0), (0, 4, 0), (0, 0, 4)])
spos_ac = np.array([[0.252, 0.15, 0.35], [0.503, 0.5, 0.5]])
s = Spline(l=0, rmax=2.0, f_g=np.array([1, 0.9, 0.1, 0.0]))
spline_aj = [[s], [s]]
c = LFC(gd, spline_aj)
c.set_positions(spos_ac)
c_ai = c.dict(zero=True)
if 1 in c_ai:
    c_ai[1][0] = 2.0
psi = gd.zeros()
c.add(psi, c_ai)

d_avv = dict([(a, np.zeros((3, 3))) for a in c.my_atom_indices])
c.second_derivative(psi, d_avv)

if 0 in d_avv:
    print(d_avv[0])

eps = 0.000001
d_aiv = c.dict(derivative=True)
pos_av = np.dot(spos_ac, gd.cell_cv)
for v in range(3):
    pos_av[0, v] += eps
    c.set_positions(np.dot(pos_av, gd.icell_cv.T))
    c.derivative(psi, d_aiv)
    if 0 in d_aiv:
        d0_v = d_aiv[0][0].copy()
    pos_av[0, v] -= 2 * eps
    c.set_positions(np.dot(pos_av, gd.icell_cv.T))
    c.derivative(psi, d_aiv)
    if 0 in d_aiv:
        d0_v -= d_aiv[0][0]
        d0_v /= -2 * eps
        print(d0_v)
        d_avv[0][v] -= d0_v
    pos_av[0, v] += eps
if 0 in d_avv:
    assert np.abs(d_avv[0]).max() < 1e-10

 

    
