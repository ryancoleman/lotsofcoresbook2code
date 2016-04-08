from __future__ import print_function
from math import sqrt, pi
import numpy as np
from gpaw.setup import create_setup
from gpaw.grid_descriptor import GridDescriptor
from gpaw.localized_functions import create_localized_functions
from gpaw.xc import XC

n = 60#40 /8 * 10
a = 10.0
gd = GridDescriptor((n, n, n), (a, a, a))
c_LL = np.identity(9, float)
a_Lg = gd.zeros(9)
nspins = 2
xc = XC('LDA')
for soft in [False]:
    s = create_setup('Cu', xc, lmax=2)
    ghat_l = s.ghat_l
    ghat_Lg = create_localized_functions(ghat_l, gd, (0.54321, 0.5432, 0.543))
    a_Lg[:] = 0.0
    ghat_Lg.add(a_Lg, c_LL)
    for l in range(3):
        for m in range(2 * l + 1):
            L = l**2 + m
            a_g = a_Lg[L]
            Q0 = gd.integrate(a_g) / sqrt(4 * pi)
            Q1_m = -gd.calculate_dipole_moment(a_g) / sqrt(4 * pi / 3)
            print(Q0)
            if l == 0:
                Q0 -= 1.0
                Q1_m[:] = 0.0
            elif l == 1:
                Q1_m[(m + 1) % 3] -= 1.0
            print(soft, l, m, Q0, Q1_m)
            assert abs(Q0) < 2e-6
            assert np.alltrue(abs(Q1_m) < 3e-5)
    b_Lg = np.reshape(a_Lg, (9, -1))
    S_LL = np.inner(b_Lg, b_Lg)
    gd.comm.sum(S_LL)
    S_LL.ravel()[::10] = 0.0
    print(max(abs(S_LL).ravel()))
    assert max(abs(S_LL).ravel()) < 3e-4
