from __future__ import print_function
from gpaw.atom.aeatom import AllElectronAtom, c
from gpaw.test import equal

Z = 79  # gold atom
kwargs = dict(ngpts=5000, alpha2=1000 * Z**2, ngauss=200)

# Test Schroedinger equation:
aea = AllElectronAtom(Z, log=None)
aea.initialize(**kwargs)

errors = []
for channel in aea.channels:
    channel.solve(-Z)
    for n in range(7):
        e = channel.e_n[n]
        e0 = -0.5 * Z**2 / (n + channel.l + 1)**2
        errors.append(abs(e / e0 - 1))
print(max(errors))
equal(max(errors), 0, 2.0e-5)

# Test Dirac equation:
aea = AllElectronAtom(Z, dirac=True, log=None)
aea.initialize(**kwargs)

errors = []
for channel in aea.channels:
    channel.solve(-Z)
    for n in range(7):
        e = channel.e_n[n]
        if channel.k > 0:
            n += 1
        e0 = (1 +
              (Z / c)**2 /
              ((channel.k**2 - (Z / c)**2)**0.5 + n)**2)**-0.5 - 1
        e0 *= c**2
        errors.append(abs(e / e0 - 1))
print(max(errors))
equal(max(errors), 0, 4.0e-5)
