# The parallel code is from Carsten Rostgaard # 

import numpy as np
from gpaw.mpi import serial_comm
from gpaw.mpi import rank, size, world
from gpaw.io import open

def set_communicator(world, rank, size, kcommsize=None):
    """Communicator inilialized."""
    # wcomm is always set to world

    wcomm = world
    
    if kcommsize is None or kcommsize == size or size == 1:
        # By default, only use parallization in kpoints
        # then kcomm is set to world communicator
        # and wS is not parallelized
        kcomm = world
        wScomm = serial_comm
        
    else:
        # If use wS parallization for storage of spectral function
        # then new kpoint and wS communicator are generated
        assert kcommsize != size
        r0 = (rank // kcommsize) * kcommsize
        ranks = np.arange(r0, r0 + kcommsize)
        kcomm = world.new_communicator(ranks)

        # wS comm generated
        r0 = rank % kcommsize
        ranks = np.arange(r0, r0+size, kcommsize)
        wScomm = world.new_communicator(ranks)

    return kcomm, wScomm, wcomm


def parallel_partition(N, commrank, commsize, reshape=True, positive=False):
    
    if reshape is True:
        if N % commsize != 0:
            N -= N % commsize
            if positive:
                N += commsize
        assert N % commsize == 0

    N_local = N // commsize
    N_residual = N - N_local * commsize
    if commrank < N_residual:
        N_local += 1
        N_start = commrank * N_local
        N_end = (commrank + 1) * N_local
    else:
        offset =  N_residual * (N_local + 1)
        N_start = offset + (commrank - N_residual) * N_local
        N_end = offset + (commrank + 1 - N_residual) * N_local
    
    return N, N_local, N_start, N_end

def parallel_partition_list(N, commrank, commsize):

    Nlist = []
    for i in range(N):
        if commrank == i % commsize:
            Nlist.append(i)

    N_local = len(Nlist)
    
    return N, N_local, Nlist

def collect_orbitals(a_xo, coords, comm, root=0):
    """Collect array distributed over orbitals to root-CPU.

    Input matrix has last axis distributed amongst CPUs,
    return is None on slaves, and the collected array on root.

    The distribution can be uneven amongst CPUs. The list coords gives the
    number of values for each CPU.
    """
    a_xo = np.ascontiguousarray(a_xo)
    if comm.size == 1:
        return a_xo

    # All slaves send their piece to ``root``:
    # There can be several sends before the corresponding receives
    # are posted, so use syncronous send here
    if comm.rank != root:
        comm.ssend(a_xo, root, 112)
        return None

    # On root, put the subdomains from the slaves into the big array
    # for the whole domain on root:
    xshape = a_xo.shape[:-1]
    Norb2 = sum(coords) # total number of orbital indices
    a_xO = np.empty(xshape + (Norb2,), a_xo.dtype)
    o = 0
    for rank, norb in enumerate(coords):
        if rank != root:
            tmp_xo = np.empty(xshape + (norb,), a_xo.dtype)
            comm.receive(tmp_xo, rank, 112)
            a_xO[..., o:o + norb] = tmp_xo
        else:
            a_xO[..., o:o + norb] = a_xo
        o += norb
    return a_xO


def collect_energies(a_ex, comm, root=0):
    """Collect array distributed over energies to root-CPU.

    Input matrix has first axis distributed evenly amongst CPUs,
    return is None on slaves, and the collected array on root.

    As the distribution is even amongst CPUs, gather can be used here.
    """
    a_ex = np.ascontiguousarray(a_ex)
    if comm.size == 1:
        return a_ex

    nenergies, xshape = a_ex.shape[0], a_ex.shape[1:]
    if comm.rank == root:
        a_Ex = np.empty((comm.size, nenergies) + xshape, a_ex.dtype)
        comm.gather(a_ex, root, a_Ex)
        return a_Ex.reshape((comm.size * nenergies,) + xshape)
    else:
        comm.gather(a_ex, root)
        return None


def SliceAlongFrequency(a_eO, coords, comm):
    """Slice along frequency axis.

    Input array has subset of energies, but full orbital matrix.
    Output has full energy, but subset of flattened orbital matrix.

    coords is a list of the number of orbital indices per cpu.
    """
    # flatten orbital axis
    a_eO = np.ascontiguousarray(a_eO).reshape(len(a_eO), -1)
    if comm.size == 1:
        return a_eO

    o = 0
    for rank, norb in enumerate(coords):
        a_eo = a_eO[:, o:o + norb].copy()
        tmp = collect_energies(a_eo, comm, root=rank)
        if rank == comm.rank:
            a_Eo = tmp
        o += norb

    return a_Eo


def SliceAlongOrbitals(a_Eo, coords, comm):
    """Slice along orbital axis.

    Input has full energy, but subset of flattened orbital matrix.
    Output array has subset of energies, but full orbital matrix.

    coords is a list of the number of orbital indices per cpu.
    """
    Norb = np.sqrt(sum(coords)).astype(int)
    nenergies = len(a_Eo) // comm.size # energy points per processor.
    a_Eo = np.ascontiguousarray(a_Eo)

    if comm.size == 1:
        return a_Eo.reshape(nenergies, Norb, Norb)

    for rank in range(comm.size):
        a_eo = a_Eo[rank * nenergies:(rank + 1) * nenergies, :].copy()
        tmp = collect_orbitals(a_eo, coords, comm, root=rank)
        if rank == comm.rank:
            a_eO = tmp.reshape(nenergies, Norb, Norb)
    return a_eO


def GatherOrbitals(a_Eo, coords, comm):
    """Slice along orbital axis.

    Input has subset of flattened orbital matrix.
    Output array has full orbital matrix.

    coords is a list of the number of orbital indices per cpu.
    """
    Norb = np.sqrt(sum(coords)).astype(int)
    if len(np.shape(a_Eo)) == 1:
        nenergies = 1
    else:
        nenergies = len(a_Eo)
    a_Eo = np.ascontiguousarray(a_Eo)

    if comm.size == 1:
        if nenergies == 1:
            return a_Eo.reshape(nenergies, Norb, Norb)[0]
        else:
            return a_Eo.reshape(nenergies, Norb, Norb)

    for rank in range(comm.size):
        tmp = collect_orbitals(a_Eo, coords, comm, root=rank)
        if rank == comm.rank:
            if nenergies==1:
                a_EO = tmp.reshape(Norb, Norb)
            else:
                a_EO = tmp.reshape(nenergies, Norb, Norb)
    return a_EO


def par_write(filename, name, comm, chi0_wGG):

    ## support only world communicator at the moment
    
    assert comm.size == size
    assert comm.rank == rank

    Nw_local, npw, npw1 = chi0_wGG.shape
    assert npw == npw1
    Nw = Nw_local * size
    
    w = open(filename, 'w', comm)
    w.dimension('Nw', Nw)
    w.dimension('npw', npw)
    w.add(name, ('Nw', 'npw', 'npw'), dtype=complex)
    if rank == 0:
        tmp = np.zeros_like(chi0_wGG[0])

    for iw in range(Nw):
        irank = iw // Nw_local
        if irank == 0:
            if rank == 0:
                w.fill(chi0_wGG[iw])
        else:
            if rank == irank:
                world.send(chi0_wGG[iw-rank*Nw_local], 0, irank+100)
            if rank == 0:
                world.receive(tmp, irank, irank+100)
                w.fill(tmp)
    if rank == 0:
        w.close()
    world.barrier()
        
def par_read(filename, name, Nw=None):

    r = open(filename, 'r')
    if Nw is None:
        Nw = r.dimension('Nw')
    else:
        assert Nw <= r.dimension('Nw')
    npw = r.dimension('npw')
    Nw_local = Nw // size
    chi0_wGG = np.zeros((Nw_local, npw, npw), dtype=complex)

    for iw in range(Nw_local):
        chi0_wGG[iw] = r.get(name, iw+rank*Nw_local)

    r.close()
    
    return chi0_wGG

def gatherv(m, N=None):

    if world.size == 1:
        return m

    ndim = m.ndim

    if  ndim == 2:
        n, N = m.shape
        assert n < N 
        M  = np.zeros((N, N), dtype=complex)
    elif ndim == 1:
        n = m.shape[0]
        M = np.zeros(N, dtype=complex)
    else:
        print('Not Implemented')
        XX
        
    n_index = np.zeros(size, dtype=int)
    world.all_gather(np.array([n]), n_index)

    root = 0
    if rank != root:
        world.ssend(m, root, 112+rank)
    else:
        for irank, n in enumerate(n_index):
            if irank == root:
                if ndim == 2:
                    M[:n_index[0] :] = m
                else:
                    M[:n_index[0]] = m
            else:
                n_start = n_index[0:irank].sum()
                n_end = n_index[0:irank+1].sum()
                if ndim == 2:
                    tmp_nN = np.zeros((n, N), dtype=complex)
                    world.receive(tmp_nN, irank, 112+irank)
                    M[n_start:n_end, :] = tmp_nN
                else:
                    tmp_n = np.zeros(n, dtype=complex)
                    world.receive(tmp_n, irank, 112+irank)
                    M[n_start:n_end] = tmp_n
    world.broadcast(M, root)
    
    return M
        
