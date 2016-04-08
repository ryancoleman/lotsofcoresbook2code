from __future__ import print_function
import numpy as np
from ase.units import Hartree
from gpaw.occupations import FermiDirac, MethfesselPaxton


class KPoint:
    eps_n = np.empty(1)
    f_n = np.empty(1)
    weight = 0.2
    s = 0
    
k = KPoint()


def f(occ, x):
    k.eps_n[0] = x
    n, dnde, x, S = occ.distribution(k, 0.0)
    return n, dnde, S

    
def test(occ):
    print(occ)
    for e in [-0.3 / Hartree, 0, 0.1 / Hartree, 1.2 / Hartree]:
        n0, d0, S0 = f(occ, e)
        x = 0.000001
        np, dp, Sp = f(occ, e + x)
        nm, dm, Sm = f(occ, e - x)
        d = -(np - nm) / (2 * x)
        dS = Sp - Sm
        dn = np - nm
        print(d - d0, dS - e * dn)
        assert abs(d - d0) < 3e-5
        assert abs(dS - e * dn) < 1e-13

for w in [0.1, 0.5]:
    test(FermiDirac(w))
    for n in range(4):
        test(MethfesselPaxton(w, n))
