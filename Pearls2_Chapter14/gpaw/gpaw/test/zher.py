from __future__ import print_function
from gpaw.utilities.blas import czher, axpy
import numpy as np
from time import time

alpha = 0.5
x = np.random.rand(3) + 1j * np.random.rand(3)
a = np.random.rand(9).reshape(3,3) + np.random.rand(9).reshape(3,3) * 1j

# make a hermitian
for i in range(3):
    for j in range(3):
        a[i,j] = a[j,i].conj()
    a[i,i] = np.real(a[i,i])

b = alpha * np.outer(x.conj(), x) + a
czher(alpha, x, a)

for i in range(3):
    for j in range(i,3):
        a[j,i] = a[i,j].conj()
        
assert np.abs(b-a).sum() < 1e-14

# testing speed
t_czher = 0
t_axpy = 0

for i in np.arange(1000):
    t0 = time()
    czher(alpha, x, a)
    t_czher += time() - t0

    t0 = time()
    xx = np.outer(x.conj(), x)
    axpy(alpha, xx, a)
    t_axpy += time() - t0

print('t_czher:', t_czher)
print('t_axpy:', t_axpy)
