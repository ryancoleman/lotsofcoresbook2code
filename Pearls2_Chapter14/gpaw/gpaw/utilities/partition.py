import numpy as np


class AtomicMatrixDistributor:
    """Class to distribute atomic dictionaries like dH_asp and D_asp."""
    def __init__(self, atom_partition, setups, kptband_comm, ns):
        self.atom_partition = atom_partition
        self.setups = setups
        self.kptband_comm = kptband_comm
        self.ns = ns
        self.new_atom_partition = self.atom_partition.to_parent_comm()

    def get_empty(self, a):
        return np.empty(self.get_shape(a))

    def get_shape(self, a):
        ni = self.setups[a].ni
        return (self.ns, ni * (ni + 1) // 2)

    def distribute(self, D_asp):
        # Right now the D are duplicated across the band/kpt comms.
        # Here we pick out a set of unique D.  With duplicates out,
        # this will be a one-to-one redistribution.
        Ddist_asp = {}
        for a in self.new_atom_partition.my_indices:
            assert self.kptband_comm.rank == 0
            Ddist_asp[a] = D_asp[a]
        self.new_atom_partition.to_even_distribution(Ddist_asp, self.get_empty)
        return Ddist_asp

    def collect(self, dHdist_asp):
        # First receive one-to-one from everywhere.
        self.new_atom_partition.from_even_distribution(dHdist_asp,
                                                       self.get_empty)
        
        # We need to broadcast across band/kpt comms now.
        # Instead of doing a broadcast for each atom, we will make a big
        # buffer for all data and broadcast that.  We just need to allocate
        # arrays of the right size.
        shapes = [self.get_shape(a) for a in self.atom_partition.my_indices]
        sizes = [np.prod(sh) for sh in shapes]
        csizes = np.cumsum(sizes)

        if len(csizes) > 0:
            bigbuf = np.empty(csizes[-1])
        else:
            bigbuf = np.empty(0) # variable name not so descriptive
        
        if self.kptband_comm.rank == 0:
            i1 = 0
            for i2, a in zip(csizes, self.atom_partition.my_indices):
                bigbuf[i1:i2] = dHdist_asp[a].ravel()
                i1 = i2
        
        self.kptband_comm.broadcast(bigbuf, 0)

        # Copy from bigbuf to reconstruct dictionary
        dH_asp = {}
        i1 = 0
        for i2, a in zip(csizes, self.atom_partition.my_indices):
            buf = self.get_empty(a)
            buf.ravel()[:] = bigbuf[i1:i2]
            dH_asp[a] = buf
            i1 = i2
        
        return dH_asp


class EvenPartitioning:
    """Represents an even partitioning of N elements over a communicator.

    For example N=17 and comm.size=5 will result in this distribution:

     * rank 0 has 3 local elements: 0, 1, 2
     * rank 1 has 3 local elements: 3, 4, 5
     * rank 2 has 3 local elements: 6, 7, 8
     * rank 3 has 4 local elements: 9, 10, 11, 12
     * rank 4 has 4 local elements: 13, 14, 15, 16

    This class uses only the 'rank' and 'size' communicator attributes."""
    def __init__(self, comm, N):
        # Conventions:
        #  n, N: local/global size
        #  i, I: local/global index
        self.comm = comm
        self.N = N
        self.nlong = -(-N // comm.size) # size of a 'long' slice
        self.nshort = N // comm.size # size of a 'short' slice
        self.longcount = N % comm.size # number of ranks with a 'long' slice
        self.shortcount = comm.size - self.longcount # ranks with 'short' slice

    def nlocal(self, rank=None):
        """Get the number of locally stored elements."""
        if rank is None:
            rank = self.comm.rank
        if rank < self.shortcount:
            return self.nshort
        else:
            return self.nlong

    def minmax(self, rank=None):
        """Get the minimum and maximum index of elements stored locally."""
        if rank is None:
            rank = self.comm.rank
        I1 = self.nshort * rank
        if rank < self.shortcount:
            I2 = I1 + self.nshort
        else:
            I1 += rank - self.shortcount
            I2 = I1 + self.nlong
        return I1, I2

    def slice(self, rank=None):
        """Get the list of indices of locally stored elements."""
        I1, I2 = self.minmax(rank=rank)
        return np.arange(I1, I2)

    def global2local(self, I):
        """Get a tuple (rank, local index) from global index I."""
        nIshort = self.nshort * self.shortcount
        if I < nIshort:
            return I // self.nshort, I % self.nshort
        else:
            Ioffset = I - nIshort
            return self.shortcount + Ioffset // self.nlong, Ioffset % self.nlong

    def local2global(self, i, rank=None):
        """Get global index I corresponding to local index i on rank."""
        if rank is None:
            rank = self.comm.rank
        return rank * self.nshort + max(rank - self.shortcount, 0) + i

    def as_atom_partition(self):
        rank_a = [self.global2local(i)[0] for i in range(self.N)]
        return AtomPartition(self.comm, rank_a)

    def get_description(self):
        lines = []
        for a in range(self.comm.size):
            elements = ', '.join(map(str, self.slice(a)))
            line = 'rank %d has %d local elements: %s' % (a, self.nlocal(a),
                                                          elements)
            lines.append(line)
        return '\n'.join(lines)


# Interface for things that can be redistributed with general_redistribute
class Redistributable:
    def get_recvbuffer(self, a): raise NotImplementedError
    def get_sendbuffer(self, a): raise NotImplementedError
    def assign(self, a): raise NotImplementedError


# Let's keep this as an independent function for now in case we change the
# classes 5 times, like we do
def general_redistribute(comm, src_rank_a, dst_rank_a, redistributable):
    # To do: it should be possible to specify duplication to several ranks
    # But how is this done best?
    requests = []
    flags = (src_rank_a != dst_rank_a)
    my_incoming_atom_indices = np.argwhere(np.bitwise_and(flags, \
        dst_rank_a == comm.rank)).ravel()
    my_outgoing_atom_indices = np.argwhere(np.bitwise_and(flags, \
        src_rank_a == comm.rank)).ravel()

    for a in my_incoming_atom_indices:
        # Get matrix from old domain:
        buf = redistributable.get_recvbuffer(a)
        requests.append(comm.receive(buf, src_rank_a[a], tag=a, block=False))
        # These arrays are not supposed to pointers into a larger,
        # contiguous buffer, so we should make a copy - except we
        # must wait until we have completed the send/receiving
        # into them, so we will do it a few lines down.
        redistributable.assign(a, buf)

    for a in my_outgoing_atom_indices:
        # Send matrix to new domain:
        buf = redistributable.get_sendbuffer(a)
        requests.append(comm.send(buf, dst_rank_a[a], tag=a, block=False))
                                  
    comm.waitall(requests)


class AtomPartition:
    """Represents atoms distributed on a standard grid descriptor."""
    def __init__(self, comm, rank_a):
        self.comm = comm
        self.rank_a = np.array(rank_a)
        self.my_indices = self.get_indices(comm.rank)
        self.natoms = len(rank_a)
    
    def reorder(self, args_a):
        # XXX use to get better load balance
        # after creating from EvenPartition
        return AtomPartition(self.comm, self.rank_a[args_a])

    def get_indices(self, rank):
        return np.where(self.rank_a == rank)[0]

    def to_parent_comm(self):
        # XXX assume communicator is strided, i.e. regular.
        # This actually imposes implicit limitations on things, but is not
        # "likely" to cause trouble with the usual communicators, i.e.
        # for gd/kd/bd.
        parent = self.comm.parent
        if parent is None:
            # This should not ordinarily be necessary, but when running with
            # AtomPAW, it is.  So let's stay out of trouble.
            return self
        
        members = self.comm.get_members()
        parent_rank_a = members[self.rank_a]
        
        # XXXX we hope and pray that our communicator is "equivalent" to
        # that which includes parent's rank0.
        assert min(members) == members[0]
        parent_rank_a -= members[0] # yuckkk
        return AtomPartition(self.comm.parent, parent_rank_a)

    def to_even_distribution(self, atomdict_ax, get_empty, copy=False):
        if copy:
            atomdict1_ax = {}
            for a, arr_x in atomdict_ax.items():
                atomdict1_ax[a] = arr_x.copy() # list-style ones??
            atomdict_ax = atomdict1_ax
        
        even_part = EvenPartitioning(self.comm,
                                     len(self.rank_a)).as_atom_partition()
        self.redistribute(even_part, atomdict_ax, get_empty)
        return atomdict_ax # XXX copy or not???

    def from_even_distribution(self, atomdict_ax, get_empty):
        even_part = EvenPartitioning(self.comm,
                                     len(self.rank_a)).as_atom_partition()
        even_part.redistribute(self, atomdict_ax, get_empty)

    def redistribute(self, new_partition, atomdict_ax, get_empty):
        assert self.comm == new_partition.comm
        # atomdict_ax may be a dictionary or a list of dictionaries

        has_many = not hasattr(atomdict_ax, 'items')
        if has_many:
            class Redist:
                def get_recvbuffer(self, a):
                    return get_empty(a)
                def assign(self, a, b_x):
                    for u, d_ax in enumerate(atomdict_ax):
                        assert a not in d_ax
                        atomdict_ax[u][a] = b_x[u]
                def get_sendbuffer(self, a):
                    return np.array([d_ax.pop(a) for d_ax in atomdict_ax])
        else:
            class Redist:
                def get_recvbuffer(self, a):
                    return get_empty(a)
                def assign(self, a, b_x):
                    assert a not in atomdict_ax
                    atomdict_ax[a] = b_x
                def get_sendbuffer(self, a):
                    return atomdict_ax.pop(a)

        general_redistribute(self.comm, self.rank_a,
                             new_partition.rank_a, Redist())
