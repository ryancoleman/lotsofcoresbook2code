"""Test of PBLAS Level 2 & 3 : rk, r2k, gemv, gemm.

The test generates random matrices A0, B0, X0, etc. on a
1-by-1 BLACS grid. They are redistributed to a mprocs-by-nprocs
BLACS grid, BLAS operations are performed in parallel, and 
results are compared against BLAS.
"""
from __future__ import print_function
import sys

import numpy as np

from gpaw.mpi import world, rank
from gpaw.test import equal
from gpaw.blacs import BlacsGrid, Redistributor, parallelprint
from gpaw.utilities import compiled_with_sl
from gpaw.utilities.blas import gemm, gemv, r2k, rk
from gpaw.utilities.scalapack import pblas_simple_gemm, pblas_simple_gemv, \
    pblas_simple_r2k, pblas_simple_rk, pblas_simple_hemm

import _gpaw

tol = 4.0e-13 # may need to be be increased if the mprocs-by-nprocs \
    # BLACS grid becomes larger

def main(M=160, N=120, K=140, seed=42, mprocs=2, nprocs=2, dtype=float):
    gen = np.random.RandomState(seed)
    grid = BlacsGrid(world, mprocs, nprocs)
    
    if (dtype==complex):
        epsilon = 1.0j
    else:
        epsilon = 0.0

    # Create descriptors for matrices on master:
    globA = grid.new_descriptor(M, K, M, K)
    globB = grid.new_descriptor(K, N, K, N)
    globC = grid.new_descriptor(M, N, M, N)
    globZ = grid.new_descriptor(K, K, K, K)
    globX = grid.new_descriptor(K, 1, K, 1)
    globY = grid.new_descriptor(M, 1, M, 1)
    globD = grid.new_descriptor(M, K, M, K)
    globS = grid.new_descriptor(M, M, M, M)
    globU = grid.new_descriptor(M, M, M, M)

    globHEC = grid.new_descriptor(K,K, K, K)
    
    # print globA.asarray()
    # Populate matrices local to master:
    A0 = gen.rand(*globA.shape) + epsilon * gen.rand(*globA.shape)
    B0 = gen.rand(*globB.shape) + epsilon * gen.rand(*globB.shape)
    D0 = gen.rand(*globD.shape) + epsilon * gen.rand(*globD.shape)
    X0 = gen.rand(*globX.shape) + epsilon * gen.rand(*globX.shape)
    
    # HEC = HEA * B
    HEA0 = gen.rand(*globHEC.shape) + epsilon * gen.rand(*globHEC.shape)
    if world.rank == 0:
        HEA0 = HEA0 + HEA0.T.conjugate() # Make H0 hermitean
    
    # Local result matrices
    Y0 = globY.empty(dtype=dtype)
    C0 = globC.zeros(dtype=dtype)
    Z0 = globZ.zeros(dtype=dtype)
    S0 = globS.zeros(dtype=dtype) # zeros needed for rank-updates
    U0 = globU.zeros(dtype=dtype) # zeros needed for rank-updates
    HEC0 = globB.zeros(dtype=dtype)

    # Local reference matrix product:
    if rank == 0:
        # C0[:] = np.dot(A0, B0)
        gemm(1.0, B0, A0, 0.0, C0)
        #gemm(1.0, A0, A0, 0.0, Z0, transa='t')
        print(A0.shape, Z0.shape)
        Z0[:] = np.dot(A0.T, A0)
        # Y0[:] = np.dot(A0, X0)
        gemv(1.0, A0, X0.ravel(), 0.0, Y0.ravel())
        r2k(1.0, A0, D0, 0.0, S0)
        rk(1.0, A0, 0.0, U0)

        HEC0[:] = np.dot(HEA0, B0)
        sM, sN = HEA0.shape
        # We don't use upper diagonal
        for i in range(sM):
           for j in range(sN):
               if i<j:
                   HEA0[i][j] = 99999.0
        if world.rank == 0:
            print(HEA0)
    assert globA.check(A0) and globB.check(B0) and globC.check(C0)
    assert globX.check(X0) and globY.check(Y0)
    assert globD.check(D0) and globS.check(S0) and globU.check(U0)

    # Create distributed destriptors with various block sizes:
    distA = grid.new_descriptor(M, K, 2, 2)
    distB = grid.new_descriptor(K, N, 2, 4)
    distC = grid.new_descriptor(M, N, 3, 2)
    distZ = grid.new_descriptor(K, K, 5, 7)
    distX = grid.new_descriptor(K, 1, 4, 1)
    distY = grid.new_descriptor(M, 1, 3, 1)
    distD = grid.new_descriptor(M, K, 2, 3)
    distS = grid.new_descriptor(M, M, 2, 2)
    distU = grid.new_descriptor(M, M, 2, 2)
    distHE = grid.new_descriptor(K, K, 2, 4)

    # Distributed matrices:
    A = distA.empty(dtype=dtype)
    B = distB.empty(dtype=dtype)
    C = distC.empty(dtype=dtype)
    Z = distZ.empty(dtype=dtype)
    X = distX.empty(dtype=dtype)
    Y = distY.empty(dtype=dtype)
    D = distD.empty(dtype=dtype)
    S = distS.zeros(dtype=dtype) # zeros needed for rank-updates
    U = distU.zeros(dtype=dtype) # zeros needed for rank-updates
    HEC = distB.zeros(dtype=dtype)
    HEA = distHE.zeros(dtype=dtype)
    Redistributor(world, globA, distA).redistribute(A0, A)
    Redistributor(world, globB, distB).redistribute(B0, B)
    Redistributor(world, globX, distX).redistribute(X0, X)
    Redistributor(world, globD, distD).redistribute(D0, D)
    Redistributor(world, globHEC, distHE).redistribute(HEA0, HEA)

    pblas_simple_gemm(distA, distB, distC, A, B, C)
    pblas_simple_gemm(distA, distA, distZ, A, A, Z, transa='T')
    pblas_simple_gemv(distA, distX, distY, A, X, Y)
    pblas_simple_r2k(distA, distD, distS, A, D, S)
    pblas_simple_rk(distA, distU, A, U)
    pblas_simple_hemm(distHE, distB, distB, HEA, B, HEC, uplo='L', side='L') 

    # Collect result back on master
    C1 = globC.empty(dtype=dtype)
    Y1 = globY.empty(dtype=dtype)
    S1 = globS.zeros(dtype=dtype) # zeros needed for rank-updates
    U1 = globU.zeros(dtype=dtype) # zeros needed for rank-updates
    HEC1 = globB.zeros(dtype=dtype)
    Redistributor(world, distC, globC).redistribute(C, C1)
    Redistributor(world, distY, globY).redistribute(Y, Y1)
    Redistributor(world, distS, globS).redistribute(S, S1)
    Redistributor(world, distU, globU).redistribute(U, U1)
    Redistributor(world, distB, globB).redistribute(HEC, HEC1)

    if rank == 0:
        gemm_err = abs(C1 - C0).max()
        gemv_err = abs(Y1 - Y0).max()
        r2k_err  = abs(S1 - S0).max()
        rk_err   = abs(U1 - U0).max()
        hemm_err = abs(HEC1 - HEC0).max()
        print('gemm err', gemm_err)
        print('gemv err', gemv_err)
        print('r2k err' , r2k_err)
        print('rk_err'  , rk_err)
        print('hemm_err', hemm_err)
    else:
        gemm_err = 0.0
        gemv_err = 0.0
        r2k_err  = 0.0
        rk_err   = 0.0
        hemm_err = 0.0

    gemm_err = world.sum(gemm_err) # We don't like exceptions on only one cpu
    gemv_err = world.sum(gemv_err)
    r2k_err  = world.sum(r2k_err)
    rk_err   = world.sum(rk_err)
    hemm_err = world.sum(hemm_err)

    equal(gemm_err, 0, tol)
    equal(gemv_err, 0, tol)
    equal(r2k_err, 0, tol)
    equal(rk_err,0, tol)
    equal(hemm_err, 0, tol)

if __name__ in ['__main__', '__builtin__']:
    if not compiled_with_sl():
        print('Not built with ScaLAPACK. Test does not apply.')
    else:
        main(dtype=float)
        main(dtype=complex)
