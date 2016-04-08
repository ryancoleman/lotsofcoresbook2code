import numpy as np
from gpaw.mpi import world, MASTER
from gpaw.blacs import BlacsGrid
from gpaw.blacs import Redistributor


def parallel_eigh(matrixfile, blacsgrid=(4, 2), blocksize=64):
    """Diagonalize matrix in parallel"""
    assert np.prod(blacsgrid) == world.size
    grid = BlacsGrid(world, *blacsgrid)

    if world.rank == MASTER:
        H_MM = np.load(matrixfile)
        assert H_MM.ndim == 2
        assert H_MM.shape[0] == H_MM.shape[1]
        NM = len(H_MM)
    else:
        NM = 0
    NM = world.sum(NM) # Distribute matrix shape to all nodes

    # descriptor for the individual blocks
    block_desc = grid.new_descriptor(NM, NM, blocksize, blocksize)

    # descriptor for global array on MASTER
    local_desc = grid.new_descriptor(NM, NM, NM, NM)

    # Make some dummy array on all the slaves
    if world.rank != MASTER:
        H_MM = local_desc.zeros()
    assert local_desc.check(H_MM)

    # The local version of the matrix
    H_mm = block_desc.empty()

    # Distribute global array to smaller blocks
    redistributor = Redistributor(world, local_desc, block_desc)
    redistributor.redistribute(H_MM, H_mm)

    # Allocate arrays for eigenvalues and -vectors
    eps_M = np.empty(NM)
    C_mm = block_desc.empty()
    block_desc.diagonalize_ex(H_mm, C_mm, eps_M)

    # Collect eigenvectors on MASTER
    C_MM = local_desc.empty()
    redistributor2 = Redistributor(world, block_desc, local_desc)
    redistributor2.redistribute(C_mm, C_MM)

    # Return eigenvalues and -vectors on Master
    if world.rank == MASTER:
        return eps_M, C_MM
    else:
        return None, None


from gpaw.utilities import compiled_with_sl
if __name__ == '__main__' and compiled_with_sl():
    # Test script which should be run on 1, 2, 4, or 8 CPUs
    
    if world.size == 1:
        blacsgrid = (1, 1)
    elif world.size == 2:
        blacsgrid = (2, 1)
    elif world.size == 4:
        blacsgrid = (2, 2)
    elif world.size == 8:
        blacsgrid = (4, 2)
    else:
        raise RuntimeError('Please use 1, 2, 4, or 8 nodes for this test')
    
    if world.rank == MASTER:
        a = np.diag(range(1,51)).astype(float)
        a.dump('H_50x50.pckl')
        eps, U = np.linalg.eigh(a)

    eps, U = parallel_eigh('H_50x50.pckl', blacsgrid, blocksize=6)
    if world.rank == MASTER:
        assert abs(eps - range(1, 51)).sum() < 1e-5
        assert abs(U - np.identity(50)).sum() < 1e-5
