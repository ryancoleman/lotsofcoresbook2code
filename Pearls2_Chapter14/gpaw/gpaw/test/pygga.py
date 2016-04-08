from ase import Atoms
from gpaw import GPAW
import numpy as np

# spin paired H2
d = 0.75
h2 = Atoms('H2', [[0, 0, 0], [0, 0, d]])
h2.center(vacuum=2.)

e = np.array([])
f = np.array([])

for XC in ['pyPBE', 'pyRPBE']:
    e = np.array([])
    f = np.array([])
    for i in [2, 0]:
        xc = XC[i:]
        calc = GPAW(nbands=-1, xc=xc, txt=None)
        h2.set_calculator(calc)
        e = np.append(e, h2.get_potential_energy())
        f = np.append(f, h2.get_forces())
        del calc

    assert np.abs(e[0] - e[1]) < 1.e-4
    assert np.sum(np.abs(f[0] - f[1])) < 1.e-10


# spin polarized O2
d = 1.2
o2 = Atoms('O2', [[0, 0, 0], [0, 0, d]], magmoms=[1., 1.])
o2.center(vacuum=2.)

for XC in ['pyPBE', 'pyRPBE']:
    e = np.array([])
    f = np.array([])
    for i in [2, 0]:
        xc = XC[i:]
        calc = GPAW(nbands=-2, xc=xc, txt=None)
        o2.set_calculator(calc)
        e = np.append(e, o2.get_potential_energy())
        f = np.append(f, o2.get_forces())
        del calc

    assert np.abs(e[0] - e[1]) < 5.e-3
    assert np.sum(np.abs(f[0] - f[1])) < 1.e-4
