from __future__ import division, print_function
from numpy.polynomial import Polynomial
import ase.units as u
from ase.lattice import bulk
from gpaw import GPAW, PW

afcc = 3.985             # Theoretical fcc lattice parameter
a = afcc * 2**(-1 / 3)   # Assuming the same volume per atom
a = afcc * (2 / 3)**0.5  # Assuming the same nearest neighbor distance

bcc = bulk('Al', 'bcc', a=a)
bcc.calc = GPAW()

# Convergence with respect to k-points:
bcc.calc.set(mode=PW(300), txt='Al-bcc-k.txt')

for k in [4, 6, 8, 10]:
    bcc.calc.set(kpts=(k, k, k))
    print(k, bcc.get_potential_energy())

# Convergence with respect to grid spacing:
bcc.calc.set(kpts=(8, 8, 8), txt='Al-bcc-ecut.txt')

for ecut in [200, 300, 400, 500]:
    bcc.calc.set(mode=PW(ecut))
    print(ecut, bcc.get_potential_energy())

# Set parameters to reasonably converged values
E = []
A = [3.0, 3.1, 3.2, 3.3]
for a in A:
    bcc = bulk('Al', 'bcc', a=a)
    bcc.calc = GPAW(mode=PW(300),
                    kpts=(8, 8, 8),
                    txt='bulk-bcc-a%.1f.txt' % a)
    E.append(bcc.get_potential_energy())

p = Polynomial.fit(A, E, 3)
a0 = p.deriv(1).roots()[0]
B = p.deriv(2)(a0) * 2 / 9 / a0 / u.J * u.m**3 * 1e-9  # GPa
print(a0, B)

assert abs(a0 - 3.1969) < 0.001
assert abs(B - 77.5) < 0.1
