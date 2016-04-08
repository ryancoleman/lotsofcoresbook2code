from __future__ import print_function
import numpy as np

from gpaw.atom.radialgd import fsbt


N = 1024
L = 50.0
h = L / N
alpha = 5.0
r = np.linspace(0, L, N, endpoint=False)
G = np.linspace(0, np.pi / h, N // 2 + 1)
n = np.exp(-alpha * r**2)

for l in range(7):
    f = fsbt(l, n * r**l, r, G)
    f0 = (np.pi**0.5 / alpha**(l + 1.5) / 2**l * G**l / 4 *
          np.exp(-G**2 / (4 * alpha)))
    print(l, abs(f-f0).max())
    assert abs(f-f0).max() < 1e-7
