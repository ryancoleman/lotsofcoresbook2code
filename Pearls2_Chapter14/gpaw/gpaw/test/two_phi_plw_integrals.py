from __future__ import print_function
import numpy as np
from math import sqrt, pi
from gpaw.lfc import LocalizedFunctionsCollection as LFC
from gpaw.grid_descriptor import GridDescriptor, RadialGridDescriptor
from gpaw.spline import Spline
from gpaw.setup import Setup
from gpaw.gaunt import gaunt as G_LLL
from gpaw.spherical_harmonics import Y
from gpaw.response.math_func import two_phi_planewave_integrals

# Initialize s, p, d (9 in total) wave and put them on grid
rc = 2.0
a = 2.5 * rc
n = 64
lmax = 2
b = 8.0
m = (lmax + 1)**2
gd = GridDescriptor([n, n, n], [a, a, a])
r = np.linspace(0, rc, 200)
g = np.exp(-(r / rc * b)**2)
splines = [Spline(l=l, rmax=rc, f_g=g) for l in range(lmax + 1)]
c = LFC(gd, [splines])
c.set_positions([(0.5, 0.5, 0.5)])
psi = gd.zeros(m)
d0 = c.dict(m)
if 0 in d0:
    d0[0] = np.identity(m)
c.add(psi, d0)

# Calculate on 3d-grid < phi_i | e**(-ik.r) | phi_j >
R_a = np.array([a/2,a/2,a/2])
rr = gd.get_grid_point_coordinates()
for dim in range(3):
    rr[dim] -= R_a[dim]

k_G = np.array([[11.,0.2,0.1],[10., 0., 10.]])
nkpt = k_G.shape[0]

d0 = np.zeros((nkpt,m,m), dtype=complex)
for i in range(m):
    for j in range(m):
        for ik in range(nkpt):
            k = k_G[ik]
            kk = np.sqrt(np.inner(k,k))
            kr = np.inner(k,rr.T).T
            expkr = np.exp(-1j * kr)
            d0[ik, i,j] = gd.integrate(psi[i] * psi[j] * expkr)

# Calculate on 1d-grid < phi_i | e**(-ik.r) | phi_j >
rgd = RadialGridDescriptor(r, np.ones_like(r) * r[1])
g = [np.exp(-(r / rc * b)**2)*r**l for l in range(lmax + 1)]
l_j = range(lmax + 1)
d1 = two_phi_planewave_integrals(k_G, rgd=rgd, phi_jg=g,
                            phit_jg=np.zeros_like(g), l_j=l_j)

d1 = d1.reshape(nkpt, m, m)

for i in range(m):
    for j in range(m):
        for ik in range(nkpt):
            if np.abs(d0[ik,i,j] - d1[ik,i,j]) > 1e-10:
                print(i, j, d0[ik,i,j]- d1[ik,i,j])

