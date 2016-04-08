"""Test correctness of vdW-DF potential."""
from math import pi
from gpaw.grid_descriptor import GridDescriptor
import numpy as np
from gpaw.test import equal
from gpaw.xc import XC
from gpaw.mpi import world

N = 8
a = 2.0
gd = GridDescriptor((N, N, N), (a, a, a))


def paired():
    xc = XC('vdW-DF')
    n = 0.3 * np.ones((1, N, N, N))
    n += 0.01 * np.cos(np.arange(N) * 2 * pi / N)
    v = 0.0 * n
    xc.calculate(gd, n, v)
    n2 = 1.0 * n
    i = 1
    n2[0, i, i, i] += 0.00002
    x = v[0, i, i, i] * gd.dv
    E2 = xc.calculate(gd, n2, v)
    n2[0, i, i, i] -= 0.00004
    E2 -= xc.calculate(gd, n2, v)
    x2 = E2 / 0.00004
    print(i, x, x2, x - x2, x / x2)
    equal(x, x2, 2e-11)


def polarized():
    xc = XC('vdW-DF')
    n = 0.04 * np.ones((2, N, N, N))
    n[1] = 0.3
    n[0] += 0.02 * np.sin(np.arange(N) * 2 * pi / N)
    n[1] += 0.2 * np.cos(np.arange(N) * 2 * pi / N)
    v = 0.0 * n
    xc.calculate(gd, n, v)
    n2 = 1.0 * n
    i = 1
    n2[0, i, i, i] += 0.00002
    x = v[0, i, i, i] * gd.dv
    E2 = xc.calculate(gd, n2, v)
    n2[0, i, i, i] -= 0.00004
    E2 -= xc.calculate(gd, n2, v)
    x2 = E2 / 0.00004
    print(i, x, x2, x - x2, x / x2)
    equal(x, x2, 1e-10)

if world.size == 1:
    polarized()
    paired()
