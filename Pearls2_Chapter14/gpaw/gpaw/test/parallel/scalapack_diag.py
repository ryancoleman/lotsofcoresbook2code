import numpy as np
from gpaw.blacs import BlacsGrid, parallelprint
from gpaw.mpi import world, rank, size
from gpaw.utilities.lapack import diagonalize
from gpaw.utilities.scalapack import scalapack_diagonalize_dc
from gpaw.blacs import Redistributor


def scal_diagonalize(A, nodes='master'):
    # Diagonalize matrix A (size N*N) with scalapack
    # Usage: eps, B = scal_diagonalize(A)
    # eps and B and the eigenvalues and eigenvectors
    # nodes = 'master': eigenvectors only available on master node
    # nodes = 'all': eigenvectors broadcast to all nodes

    # make sure A is N*N, and hermitian
    N = A.shape[0]
    assert A.shape[0] == A.shape[1]
    for i in range(N):
        for j in range(i, N):
            assert A[i,j] == A[j,i].conj()

    # create blacs descriptor
    mb = 64
    g = BlacsGrid(world, 2, size//2)
    nndesc1 = g.new_descriptor(N, N, N,  N) 
    nndesc2 = g.new_descriptor(N, N, mb, mb)

    # distribute A to blacs grid A_
    if rank != 0:
        A = nndesc1.zeros(dtype=A.dtype)
    A_ = nndesc2.empty(dtype=A.dtype)
    redistributor = Redistributor(world, nndesc1, nndesc2)
    redistributor.redistribute(A, A_)

    # diagonalize
    B_ = nndesc2.zeros(dtype=A.dtype)
    eps = np.zeros(N,dtype=A.dtype)
    nndesc2.diagonalize_dc(A_, B_, eps, 'L')

    # distribute the eigenvectors to master
    B = np.zeros_like(A)
    redistributor = Redistributor(world, nndesc2, nndesc1)
    redistributor.redistribute(B_, B)

    if nodes == 'master':
        return eps, B
    elif nodes == 'all':
        if rank != 0:
            B = np.zeros((N, N))
        world.broadcast(B, 0)
        return eps, B
    

# generate a matrix
N = 512
A = np.arange(N**2,dtype=float).reshape(N,N)
for i in range(N):
    for j in range(i,N):
        A[i,j] = A[j,i]

# diagonalize
eps, B = scal_diagonalize(A)


check = 1
if check and rank == 0:
# check whether it gives the same result with lapack
    eps1 = np.zeros(N)
    diagonalize(A, eps1)
    assert np.abs(eps-eps1).sum() < 1e-6
    
    for i in range(N//size):
        # the eigenvectors are row of the matrix, it can be differ by a minus sign.
        if np.abs(A[i,:] - B[i,:]).sum() > 1e-6:
            if np.abs(A[i,:] + B[i,:]).sum() > 1e-6:
                raise ValueError('Check !')

