from __future__ import print_function
import numpy as np

# Set-up a simple matrix in parallel, diagonalize using ScaLAPACK
# D&C driver then compare *eigenvalues* with serial LAPACK diagonlize
from time import time

from gpaw.blacs import BlacsGrid, parallelprint
from gpaw.mpi import world, rank, size
from gpaw.utilities import compiled_with_sl
from gpaw.utilities.lapack import diagonalize
from gpaw.utilities.scalapack import scalapack_diagonalize_dc, scalapack_set, \
    scalapack_zero

switch_uplo = {'U': 'L', 'L': 'U'}

rank = world.rank

# tol
tol = 5e-12 # eigenvalue tolerance

def main(nbands=1000, mprocs=2, mb=64):
    # Set-up BlacsGrud
    grid = BlacsGrid(world, mprocs, mprocs)

    # Create descriptor
    nndesc = grid.new_descriptor(nbands, nbands, mb, mb)
    H_nn = nndesc.empty(dtype=float) # outside the BlacsGrid these are size zero
    C_nn = nndesc.empty(dtype=float) # outside the BlacsGrid these are size zero
    eps_N  = np.empty((nbands), dtype=float) # replicated array on all MPI tasks 
    # Fill ScaLAPACK array
    alpha = 0.1 # off-diagonal
    beta = 75.0 # diagonal
    uplo = 'L' # lower-triangular
    scalapack_set(nndesc, H_nn, alpha, beta, uplo)
    scalapack_zero(nndesc, H_nn, switch_uplo[uplo])

    t1 = time()
    # either interface will work, we recommend use the latter interface
    # scalapack_diagonalize_dc(nndesc, H_nn.copy(), C_nn, eps_N, 'L')
    nndesc.diagonalize_dc(H_nn.copy(), C_nn, eps_N) 
    t2 = time()
    world.broadcast(eps_N, 0) # all MPI tasks now have eps_N
    world.barrier() # wait for everyone to finish

    if rank == 0:
        print('ScaLAPACK diagonalize_dc', t2-t1)

    # Create replicated NumPy array    
    diagonal = np.eye(nbands,dtype=float)
    offdiagonal = np.tril(np.ones((nbands,nbands)), -1)
    H0 = beta*diagonal + alpha*offdiagonal
    E0 = np.empty((nbands), dtype=float)

    t1 = time()
    diagonalize(H0,E0)
    t2 = time()

    if rank == 0:
        print('LAPACK diagonalize', t2-t1)

    delta = abs(E0-eps_N).max()
    if rank == 0:
        print(delta)
        assert delta < tol

if __name__ in ['__main__', '__builtin__']:
    if not compiled_with_sl():
        print('Not built with ScaLAPACK. Test does not apply.')
    else:
        main()
