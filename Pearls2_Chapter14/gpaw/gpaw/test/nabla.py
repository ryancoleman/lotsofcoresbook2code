from __future__ import print_function
import numpy as np
from gpaw.lfc import LocalizedFunctionsCollection as LFC
from gpaw.grid_descriptor import GridDescriptor
from gpaw.atom.radialgd import EquidistantRadialGridDescriptor
from gpaw.spline import Spline
from gpaw.setup import Setup

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
c.set_positions([(0, 0, 0)])
psi = gd.zeros(m)
d0 = c.dict(m)
if 0 in d0:
    d0[0] = np.identity(m)
c.add(psi, d0)
d1 = c.dict(m, derivative=True)
c.derivative(psi, d1)
class TestSetup(Setup):
    l_j = range(lmax + 1)
    nj = lmax + 1
    ni = m
    def __init__(self):
        pass
rgd = EquidistantRadialGridDescriptor(r[1], len(r))
g = [np.exp(-(r / rc * b)**2) * r**l for l in range(lmax + 1)]
d2 = TestSetup().get_derivative_integrals(rgd, g, np.zeros_like(g))
if 0 in d1:
    print(abs(d1[0] - d2).max())
    assert abs(d1[0] - d2).max() < 2e-6



