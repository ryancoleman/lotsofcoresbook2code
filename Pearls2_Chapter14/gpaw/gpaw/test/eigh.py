from __future__ import print_function
import numpy as np
tol = 1e-8
a = np.eye(3, dtype=complex)
a[1:, 0] = 0.01j
w0 = [0.98585786, 1.0, 1.01414214]
# NumPy's Diagonalize
from numpy.linalg import eigh
w = eigh(a)[0]
print(w)
assert abs(w - w0).max() < tol
# LAPACK's QR Diagonalize
from gpaw.utilities.lapack import diagonalize
diagonalize(a.copy(), w)
print(w)
assert abs(w - w0).max() < tol
# LAPACK's MR3 Diagonalize
# Requires Netlib LAPACK 3.2.1 or later
# from gpaw.utilities.lapack import diagonalize_mr3
# z = np.zeros_like(a)
# diagonalize_mr3(a, w, z)
# print w
# assert abs(w - w0).max() < tol

