from __future__ import print_function

import numpy as np

from gpaw.mpi import world
from gpaw.blacs import BlacsGrid, Redistributor
from gpaw.utilities import compiled_with_sl

def test(comm, M, N, mcpus, ncpus, mb, nb):
    grid0 = BlacsGrid(comm, 1, 1)
    desc0 = grid0.new_descriptor(M, N, M, N, 0, 0)
    A_mn = desc0.zeros(dtype=float)
    A_mn[:] = comm.size + 1

    grid1 = BlacsGrid(comm, mcpus, ncpus)
    desc1 = grid1.new_descriptor(M, N, mb, nb, 0, 0) # ???
    B_mn = desc1.zeros(dtype=float)
    B_mn[:] = comm.rank

    if comm.rank == 0:
        msg = 'Slices of global matrix indices by rank'
        print(msg)
        print('-' * len(msg))

    for rank in range(comm.size):
        comm.barrier()
        if rank == comm.rank:
            print('Rank %d:' % rank)
            last_Mstart = -1
            for Mstart, Mstop, Nstart, Nstop, block in desc1.my_blocks(B_mn):
                if Mstart > last_Mstart and last_Mstart >= 0:
                    print()
                print('[%3d:%3d, %3d:%3d]' % (Mstart, Mstop, Nstart, Nstop), end=' ')
                last_Mstart = Mstart
                assert (block == comm.rank).all()
                #print block
                #print
            print()
            print()
        comm.barrier()
    
    redistributor = Redistributor(comm, desc1, desc0)
    redistributor.redistribute(B_mn, A_mn)

    if comm.rank == 0:
        msg = 'Rank where each element of the global matrix is stored'
        print(msg)
        print('-' * len(msg))
        print(A_mn)

if __name__ in ['__main__', '__builtin__']:
    if not compiled_with_sl():
        print('Not built with ScaLAPACK. Test does not apply.')
    else:
        M, N = 10, 10
        mb, nb = 2, 2
        mcpus = int(np.ceil(world.size**0.5))
        ncpus = world.size // mcpus

        if world.rank == 0:
            print('world size:   ', world.size)
            print('M x N:        ', M, 'x', N)
            print('mcpus x ncpus:', mcpus, 'x', ncpus)
            print('mb x nb:      ', mb, 'x', nb)
            print()

        test(world, M, N, mcpus, ncpus, mb, nb)
