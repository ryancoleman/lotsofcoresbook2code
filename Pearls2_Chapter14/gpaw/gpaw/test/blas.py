import numpy as np
from gpaw.utilities.blas import \
     gemm, axpy, r2k, rk, gemmdot, rotate, dotc, dotu
from gpaw.utilities.tools import tri2full

a = np.arange(5 * 7).reshape(5, 7) + 4.
a2 = np.arange(3 * 7).reshape(3, 7) + 3.
b = np.arange(7) - 2.

# Check gemmdot with floats
assert np.all(np.dot(a, b) == gemmdot(a, b))
assert np.all(np.dot(a, a2.T) == gemmdot(a, a2, trans='t'))
assert np.all(np.dot(a, a2.T) == gemmdot(a, a2, trans='c'))
assert np.dot(b, b) == gemmdot(b, b)

# Check gemmdot with complex arrays
a = a * (2 + 1.j)
a2 = a2 * (-1 + 3.j)
b = b * (3 - 2.j)
assert np.all(np.dot(a, b) == gemmdot(a, b))
assert np.all(np.dot(a, a2.T) == gemmdot(a, a2, trans='t'))
assert np.all(np.dot(a, a2.T.conj()) == gemmdot(a, a2, trans='c'))
assert np.dot(b, b) == gemmdot(b, b, trans='n')
assert np.dot(b, b.conj()) == gemmdot(b, b, trans='c')
assert np.vdot(a, 5.j * a) == dotc(a, 5.j * a)

# Check gemm for transa='n'
a2 = np.arange(7 * 5 * 1 * 3).reshape(7, 5, 1, 3) * (-1. + 4.j) + 3.
c = np.tensordot(a, a2, [1, 0])
gemm(1., a2, a, -1., c, 'n')
assert not c.any()

# Check gemm for transa='c'
a = np.arange(4 * 5 * 1 * 3).reshape(4, 5, 1, 3) * (3. - 2.j) + 4.
c = np.tensordot(a, a2.conj(), [[1, 2, 3], [1, 2, 3]])
gemm(1., a2, a, -1., c, 'c')
assert not c.any()

# Check axpy
c = 5.j * a
axpy(-5.j, a, c)
assert not c.any()

# Check rk
c = np.tensordot(a, a.conj(), [[1, 2, 3], [1, 2, 3]])
rk(1., a, -1., c)
tri2full(c)
assert not c.any()

# Check gemmdot for transa='c'
c = np.tensordot(a, a2.conj(), [-1, -1])
gemmdot(a, a2, beta=-1., out=c, trans='c')
assert not c.any()

# Check gemmdot for transa='n'
a2.shape = 3, 7, 5, 1
c = np.tensordot(a, a2, [-1, 0])
gemmdot(a, a2, beta=-1., out=c, trans='n')
assert not c.any()

# Check r2k
a2 = 5. * a
c = np.tensordot(a, a2.conj(), [[1, 2, 3], [1, 2, 3]])
r2k(.5, a, a2, -1., c)
tri2full(c)
assert not c.any()
