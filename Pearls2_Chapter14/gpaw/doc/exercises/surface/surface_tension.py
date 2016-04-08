from __future__ import print_function
from gpaw.io.array import load_array

N, E = load_array('e6x6.dat', transpose=True)

A = 4.05**2 / 2
Ecoh = 3.34

sigma_emt = (1 - (8 / 12.)**.5) * Ecoh / A
sigma_dft = (N[1:] * E[:-1] - N[:-1] * E[1:]) / (2 * A)

print('EMT: %4.1f meV/A**2' % (sigma_emt * 1000))
for n, s in zip(N[1:], sigma_dft):
    print(' %i %6.1f meV/A**2' % (n, s * 1000))

# EMT: 74.7 meV/A**2
#  2   77.2 meV/A**2
#  3  115.2 meV/A**2
#  4   60.1 meV/A**2
#  5   43.4 meV/A**2
#  6  101.3 meV/A**2
