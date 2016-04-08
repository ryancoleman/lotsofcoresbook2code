from math import pi
from gpaw.grid_descriptor import GridDescriptor
from gpaw.xc import XC
from gpaw.xc.noncollinear import NonCollinearLDAKernel, NonCollinearFunctional
import numpy as np
from gpaw.test import equal

N = 20
a = 1.0
gd = GridDescriptor((N, N, N), (a, a, a))

for xc in [XC('LDA'), XC('PBE')]:
    n = gd.empty(2)
    n.fill(0.03)
    z = np.arange(gd.beg_c[2], gd.end_c[2]) * a / N
    n[:] += 0.01 * np.sin(2 * pi * z / a)
    n[1] += 0.02 + 0.01 * np.cos(2 * pi * z / a)
    n /= 2

    v = 0.0 * n
    E = xc.calculate(gd, n, v)

    here = (gd.beg_c[0] <= 1 < gd.end_c[0] and
            gd.beg_c[1] <= 2 < gd.end_c[1] and
            gd.beg_c[2] <= 3 < gd.end_c[2])
    if here:
        x = v[-1, 1, 2, 3] * gd.dv
        n[-1, 1, 2, 3] += 0.000001
    Ep = xc.calculate(gd, n, v)
    if here:
        n[-1, 1, 2, 3] -= 0.000002
    Em = xc.calculate(gd, n, v)
    x2 = (Ep - Em) / 0.000002
    if here:
        print(xc.name, E, x, x2, x - x2)
        equal(x, x2, 1e-11)
        n[-1, 1, 2, 3] += 0.000001

    if 0:#xc.type == 'LDA':
        xc = XC(NonCollinearLDAKernel())
    else:
        xc = NonCollinearFunctional(xc)

    n2 = gd.zeros(4)
    n2[0] = n.sum(0)
    n2[3] = n[0] - n[1]
    E2 = xc.calculate(gd, n2)
    print(E, E2-E)
    assert abs(E2 - E) < 1e-11
    n2[1] = 0.1 * n2[3]
    n2[2] = 0.2 * n2[3]
    n2[3] *= (1 - 0.1**2 - 0.2**2)**0.5
    v = n2 * 0
    E2 = xc.calculate(gd, n2, v)
    print(E, E2-E)
    assert abs(E2 - E) < 1e-11

    for i in range(4):
        if here:
            x = v[i, 1, 2, 3] * gd.dv
            n2[i, 1, 2, 3] += 0.000001
        Ep = xc.calculate(gd, n2)
        if here:
            n2[i, 1, 2, 3] -= 0.000002
        Em = xc.calculate(gd, n2)
        x2 = (Ep - Em) / 0.000002
        if here:
            print(i, x, x2, x - x2)
            equal(x, x2, 1e-11)
            n2[i, 1, 2, 3] += 0.000001
