import numpy as np
from gpaw.lfc import LocalizedFunctionsCollection as LFC
from gpaw.grid_descriptor import GridDescriptor
from gpaw.spline import Spline
a = 4.0
gd = GridDescriptor(N_c=[16, 20, 20], cell_cv=[a, a + 1, a + 2],
                    pbc_c=(0, 1, 1))
spos_ac = np.array([[0.25, 0.15, 0.35], [0.5, 0.5, 0.5]])
kpts_kc = None
s = Spline(l=0, rmax=2.0, f_g=np.array([1, 0.9, 0.1, 0.0]))
p = Spline(l=1, rmax=2.0, f_g=np.array([1, 0.9, 0.1, 0.0]))
spline_aj = [[s], [s, p]]
c = LFC(gd, spline_aj, cut=True, forces=True)
c.set_positions(spos_ac)
C_ani = c.dict(3, zero=True)
if 1 in C_ani:
    C_ani[1][:, 1:] = np.eye(3)
psi = gd.zeros(3)
c.add(psi, C_ani)
c.integrate(psi, C_ani)
if 1 in C_ani:
    d = C_ani[1][:, 1:].diagonal()
    assert d.ptp() < 4e-6
    C_ani[1][:, 1:] -= np.diag(d)
    assert abs(C_ani[1]).max() < 5e-17
d_aniv = c.dict(3, derivative=True)
c.derivative(psi, d_aniv)
if 1 in d_aniv:
    for v in range(3):
        assert abs(d_aniv[1][v - 1, 0, v] + 0.2144) < 5e-5
        d_aniv[1][v - 1, 0, v] = 0
    assert abs(d_aniv[1]).max() < 3e-16
eps = 0.0001
pos_av = np.dot(spos_ac, gd.cell_cv)
for v in range(3):
    pos_av[0, v] += eps
    c.set_positions(np.dot(pos_av, gd.icell_cv.T))
    c.integrate(psi, C_ani)
    if 0 in d_aniv:
        C0_n = C_ani[0][:, 0].copy()
    pos_av[0, v] -= 2 * eps
    c.set_positions(np.dot(pos_av, gd.icell_cv.T))
    c.integrate(psi, C_ani)
    if 0 in d_aniv:
        C0_n -= C_ani[0][:, 0]
        C0_n /= -2 * eps
        assert abs(C0_n - d_aniv[0][:, 0, v]).max() < 1e-8
    pos_av[0, v] += eps
 

    
