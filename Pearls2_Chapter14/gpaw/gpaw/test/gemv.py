from __future__ import print_function

import time
import numpy as np
from gpaw.utilities.blas import gemmdot, dotu, gemv

def getrand(shape, dtype):
    if type(shape) is int:
        shape = (shape,)
    nelements = np.prod(shape)
    if dtype == float:
        return np.random.normal(size=nelements).reshape(shape)
    elif dtype == complex:
        return (np.random.uniform(size=nelements) * np.exp(1j \
            * np.random.uniform(0,2*np.pi,size=nelements))).reshape(shape)
    else:
        raise ValueError('Unsupported dtype "%s".' % dtype)

P = 12
Q = 1456
L = 378
G = 5013

dtype = float
beta = 0.0 # note that beta=0.0 and beta=1.0 have different gemv subroutines
refperf = 1e8 # each test should last roughly 1 s with this performance
test_gemmdot = False # gemmdot is unbelievable slow for matrix-vector products

mem = 0
itemsize = np.nbytes[np.dtype(dtype)]
mem += P*Q*L*itemsize #B_pqL
mem += L*itemsize #Y_L
mem += P*Q*itemsize #BY_pq
mem += Q*G*itemsize #n_qg
mem += G*itemsize #x_g
mem += Q*itemsize #nx_q

print('Estimated memory: %8.5f MB' % (mem/1024**2.,))

B_pqL = getrand((P,Q,L), dtype)
Y_L = getrand(L, dtype)
n_qg = getrand((Q,G), dtype)
x_g = getrand(G, dtype)

# -------------------------------------------------------------------

print('\n%s\nBY_pq calculations\n%s\n' % ('='*40, '='*40))
numflop = (dtype==float and 2 or 8)*P*Q*L
numreps = 1+int(refperf/numflop)

# Reference value
BY0_pq = np.dot(B_pqL, Y_L)

t = time.time()
for n in range(numreps):
    BY1_pq = np.dot(B_pqL, Y_L)
t = time.time()-t
performance = numflop*numreps/t
print('dot    : %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
assert np.abs(BY0_pq-BY1_pq).max()<5e-12
del BY1_pq

if test_gemmdot:
    BY2_pq = np.empty((P,Q), dtype)
    t = time.time()
    for n in range(numreps):
        BY2_pq.fill(0.0)
        gemmdot(B_pqL, Y_L, 1.0, beta, BY2_pq)
    t = time.time()-t
    performance = numflop*numreps/t
    print('gemmdot: %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
    assert np.abs(BY0_pq-BY2_pq).max()<5e-12
    del BY2_pq

BY3_pq = np.empty((P,Q), dtype)
t = time.time()
for n in range(numreps):
    BY3_pq.fill(0.0)
    gemv(1.0, B_pqL, Y_L, beta, BY3_pq, 't')
t = time.time()-t
performance = numflop*numreps/t
print('gemvT  : %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
assert np.abs(BY0_pq-BY3_pq).max()<5e-12
del BY3_pq

B_xL = B_pqL.reshape((P*Q,L))
BY4_x = np.empty(P*Q, dtype)
t = time.time()
for n in range(numreps):
    BY4_x.fill(0.0)
    gemv(1.0, B_xL, Y_L, beta, BY4_x, 't')
t = time.time()-t
performance = numflop*numreps/t
print('gemvT2D: %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
assert np.abs(BY0_pq-BY4_x.reshape((P,Q))).max()<5e-12
del B_xL, BY4_x

BT_Lqp = B_pqL.T.copy()
BY5T_qp = np.empty((Q,P), dtype)
t = time.time()
for n in range(numreps):
    BY5T_qp.fill(0.0)
    gemv(1.0, BT_Lqp, Y_L, beta, BY5T_qp, 'n')
t = time.time()-t
performance = numflop*numreps/t
print('gemvN  : %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
assert np.abs(BY0_pq-BY5T_qp.T).max()<5e-12
del BT_Lqp, BY5T_qp

BT_Lx = B_pqL.T.reshape((L,Q*P)).copy()
BY6T_x = np.empty(Q*P, dtype)
t = time.time()
for n in range(numreps):
    BY6T_x.fill(0.0)
    gemv(1.0, BT_Lx, Y_L, beta, BY6T_x, 'n')
t = time.time()-t
performance = numflop*numreps/t
print('gemvN2D: %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
assert np.abs(BY0_pq-BY6T_x.reshape((Q,P)).T).max()<5e-12
del BT_Lx, BY6T_x

# -------------------------------------------------------------------

print('\n%s\nnx_q calculations\n%s\n' % ('='*40, '='*40))
numflop = (dtype==float and 2 or 8)*Q*G
numreps = 1+int(refperf/numflop)

# Reference value
nx0_q = np.dot(n_qg, x_g)

t = time.time()
for n in range(numreps):
    nx1_q = np.dot(n_qg, x_g)
t = time.time()-t
performance = numflop*numreps/t
print('dot    : %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
assert np.abs(nx0_q-nx1_q).max()<5e-12
del nx1_q

if test_gemmdot:
    nx2_q = np.empty(Q, dtype)
    t = time.time()
    for n in range(numreps):
        nx2_q.fill(0.0)
        gemmdot(n_qg, x_g, 1.0, beta, nx2_q)
    t = time.time()-t
    performance = numflop*numreps/t
    print('gemmdot: %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
    assert np.abs(nx0_q-nx2_q).max()<5e-12
    del nx2_q

nx3_q = np.empty(Q, dtype)
t = time.time()
for n in range(numreps):
    nx3_q.fill(0.0)
    gemv(1.0, n_qg, x_g, beta, nx3_q, 't')
t = time.time()-t
performance = numflop*numreps/t
print('gemvT  : %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
assert np.abs(nx0_q-nx3_q).max()<5e-12
del nx3_q

nT_gq = n_qg.T.copy()
nx4_q = np.empty(Q, dtype)
t = time.time()
for n in range(numreps):
    nx4_q.fill(0.0)
    gemv(1.0, nT_gq, x_g, beta, nx4_q, 'n')
t = time.time()-t
performance = numflop*numreps/t
print('gemvN  : %8.5f s, %8.5f Mflops' % (t,performance/1024**2.))
assert np.abs(nx0_q-nx4_q).max()<5e-12
del nT_gq, nx4_q

