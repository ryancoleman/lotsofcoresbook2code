# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

import sys
import time
import atexit
import pickle
import numpy as np

from gpaw import debug
from gpaw import dry_run as dry_run_size
from gpaw.utilities import is_contiguous
from gpaw.utilities import gcd
from gpaw.utilities.tools import md5_array

import _gpaw

MASTER = 0


class _Communicator:
    def __init__(self, comm, parent=None):
        """Construct a wrapper of the C-object for any MPI-communicator.

        Parameters:

        comm: MPI-communicator
            Communicator.

        Attributes:

        ============  ======================================================
        ``size``      Number of ranks in the MPI group.
        ``rank``      Number of this CPU in the MPI group.
        ``parent``    Parent MPI-communicator.
        ============  ======================================================
        """
        self.comm = comm
        self.size = comm.size
        self.rank = comm.rank
        self.parent = parent  # XXX check C-object against comm.parent?

    def new_communicator(self, ranks):
        """Create a new MPI communicator for a subset of ranks in a group.
        Must be called with identical arguments by all relevant processes.

        Note that a valid communicator is only returned to the processes
        which are included in the new group; other ranks get None returned.

        Parameters:

        ranks: ndarray (type int)
            List of integers of the ranks to include in the new group.
            Note that these ranks correspond to indices in the current
            group whereas the rank attribute in the new communicators
            correspond to their respective index in the subset.

        """

        comm = self.comm.new_communicator(ranks)
        if comm is None:
            # This cpu is not in the new communicator:
            return None
        else:
            return _Communicator(comm, parent=self)

    def sum(self, a, root=-1):
        """Perform summation by MPI reduce operations of numerical data.

        Parameters:

        a: ndarray or value (type int, float or complex)
            Numerical data to sum over all ranks in the communicator group.
            If the data is a single value of type int, float or complex,
            the result is returned because the input argument is immutable.
            Otherwise, the reduce operation is carried out in-place such
            that the elements of the input array will represent the sum of
            the equivalent elements across all processes in the group.
        root: int (default -1)
            Rank of the root process, on which the outcome of the reduce
            operation is valid. A root rank of -1 signifies that the result
            will be distributed back to all processes, i.e. a broadcast.

        """
        if isinstance(a, (int, float, complex)):
            return self.comm.sum(a, root)
        else:
            tc = a.dtype
            assert tc == int or tc == float or tc == complex
            assert is_contiguous(a, tc)
            assert root == -1 or 0 <= root < self.size
            self.comm.sum(a, root)

    def product(self, a, root=-1):
        """Do multiplication by MPI reduce operations of numerical data.

        Parameters:

        a: ndarray or value (type int or float)
            Numerical data to multiply across all ranks in the communicator
            group. NB: Find the global product from the local products.
            If the data is a single value of type int or float (no complex),
            the result is returned because the input argument is immutable.
            Otherwise, the reduce operation is carried out in-place such
            that the elements of the input array will represent the product
            of the equivalent elements across all processes in the group.
        root: int (default -1)
            Rank of the root process, on which the outcome of the reduce
            operation is valid. A root rank of -1 signifies that the result
            will be distributed back to all processes, i.e. a broadcast.

        """
        if isinstance(a, (int, float)):
            return self.comm.product(a, root)
        else:
            tc = a.dtype
            assert tc == int or tc == float
            assert is_contiguous(a, tc)
            assert root == -1 or 0 <= root < self.size
            self.comm.product(a, root)

    def max(self, a, root=-1):
        """Find maximal value by an MPI reduce operation of numerical data.

        Parameters:

        a: ndarray or value (type int or float)
            Numerical data to find the maximum value of across all ranks in
            the communicator group. NB: Find global maximum from local max.
            If the data is a single value of type int or float (no complex),
            the result is returned because the input argument is immutable.
            Otherwise, the reduce operation is carried out in-place such
            that the elements of the input array will represent the max of
            the equivalent elements across all processes in the group.
        root: int (default -1)
            Rank of the root process, on which the outcome of the reduce
            operation is valid. A root rank of -1 signifies that the result
            will be distributed back to all processes, i.e. a broadcast.

        """
        if isinstance(a, (int, float)):
            return self.comm.max(a, root)
        else:
            tc = a.dtype
            assert tc == int or tc == float
            assert is_contiguous(a, tc)
            assert root == -1 or 0 <= root < self.size
            self.comm.max(a, root)

    def min(self, a, root=-1):
        """Find minimal value by an MPI reduce operation of numerical data.

        Parameters:

        a: ndarray or value (type int or float)
            Numerical data to find the minimal value of across all ranks in
            the communicator group. NB: Find global minimum from local min.
            If the data is a single value of type int or float (no complex),
            the result is returned because the input argument is immutable.
            Otherwise, the reduce operation is carried out in-place such
            that the elements of the input array will represent the min of
            the equivalent elements across all processes in the group.
        root: int (default -1)
            Rank of the root process, on which the outcome of the reduce
            operation is valid. A root rank of -1 signifies that the result
            will be distributed back to all processes, i.e. a broadcast.

        """
        if isinstance(a, (int, float)):
            return self.comm.min(a, root)
        else:
            tc = a.dtype
            assert tc == int or tc == float
            assert is_contiguous(a, tc)
            assert root == -1 or 0 <= root < self.size
            self.comm.min(a, root)

    def scatter(self, a, b, root):
        """Distribute data from one rank to all other processes in a group.

        Parameters:

        a: ndarray (ignored on all ranks different from root; use None)
            Source of the data to distribute, i.e. send buffer on root rank.
        b: ndarray
            Destination of the distributed data, i.e. local receive buffer.
            The size of this array multiplied by the number of process in
            the group must match the size of the source array on the root.
        root: int
            Rank of the root process, from which the source data originates.

        The reverse operation is ``gather``.

        Example::

          # The master has all the interesting data. Distribute it.
          if comm.rank == 0:
              data = np.random.normal(size=N*comm.size)
          else:
              data = None
          mydata = np.empty(N, dtype=float)
          comm.scatter(data, mydata, 0)

          # .. which is equivalent to ..

          if comm.rank == 0:
              # Extract my part directly
              mydata[:] = data[0:N]
              # Distribute parts to the slaves
              for rank in range(1, comm.size):
                  buf = data[rank*N:(rank+1)*N]
                  comm.send(buf, rank, tag=123)
          else:
              # Receive from the master
              comm.receive(mydata, 0, tag=123)

        """
        if self.rank == root:
            assert a.dtype == b.dtype
            assert a.size == self.size * b.size
            assert a.flags.contiguous
        assert b.flags.contiguous
        assert 0 <= root < self.size
        self.comm.scatter(a, b, root)

    def all_gather(self, a, b):
        """Gather data from all ranks onto all processes in a group.

        Parameters:

        a: ndarray
            Source of the data to gather, i.e. send buffer of this rank.
        b: ndarray
            Destination of the distributed data, i.e. receive buffer.
            The size of this array must match the size of the distributed
            source arrays multiplied by the number of process in the group.

        Example::

          # All ranks have parts of interesting data. Gather on all ranks.
          mydata = np.random.normal(size=N)
          data = np.empty(N*comm.size, dtype=float)
          comm.all_gather(mydata, data)

          # .. which is equivalent to ..

          if comm.rank == 0:
              # Insert my part directly
              data[0:N] = mydata
              # Gather parts from the slaves
              buf = np.empty(N, dtype=float)
              for rank in range(1, comm.size):
                  comm.receive(buf, rank, tag=123)
                  data[rank*N:(rank+1)*N] = buf
          else:
              # Send to the master
              comm.send(mydata, 0, tag=123)
          # Broadcast from master to all slaves
          comm.broadcast(data, 0)

        """
        assert a.flags.contiguous
        assert b.flags.contiguous
        assert b.dtype == a.dtype
        assert (b.shape[0] == self.size and a.shape == b.shape[1:] or
                a.size * self.size == b.size)
        self.comm.all_gather(a, b)

    def gather(self, a, root, b=None):
        """Gather data from all ranks onto a single process in a group.

        Parameters:

        a: ndarray
            Source of the data to gather, i.e. send buffer of this rank.
        root: int
            Rank of the root process, on which the data is to be gathered.
        b: ndarray (ignored on all ranks different from root; default None)
            Destination of the distributed data, i.e. root's receive buffer.
            The size of this array must match the size of the distributed
            source arrays multiplied by the number of process in the group.

        The reverse operation is ``scatter``.

        Example::

          # All ranks have parts of interesting data. Gather it on master.
          mydata = np.random.normal(size=N)
          if comm.rank == 0:
              data = np.empty(N*comm.size, dtype=float)
          else:
              data = None
          comm.gather(mydata, 0, data)

          # .. which is equivalent to ..

          if comm.rank == 0:
              # Extract my part directly
              data[0:N] = mydata
              # Gather parts from the slaves
              buf = np.empty(N, dtype=float)
              for rank in range(1, comm.size):
                  comm.receive(buf, rank, tag=123)
                  data[rank*N:(rank+1)*N] = buf
          else:
              # Send to the master
              comm.send(mydata, 0, tag=123)

        """
        assert a.flags.contiguous
        assert 0 <= root < self.size
        if root == self.rank:
            assert b.flags.contiguous and b.dtype == a.dtype
            assert (b.shape[0] == self.size and a.shape == b.shape[1:] or
                    a.size * self.size == b.size)
            self.comm.gather(a, root, b)
        else:
            assert b is None
            self.comm.gather(a, root)

    def broadcast(self, a, root):
        """Share data from a single process to all ranks in a group.

        Parameters:

        a: ndarray
            Data, i.e. send buffer on root rank, receive buffer elsewhere.
            Note that after the broadcast, all ranks have the same data.
        root: int
            Rank of the root process, from which the data is to be shared.

        Example::

          # All ranks have parts of interesting data. Take a given index.
          mydata[:] = np.random.normal(size=N)

          # Who has the element at global index 13? Everybody needs it!
          index = 13
          root, myindex = divmod(index, N)
          element = np.empty(1, dtype=float)
          if comm.rank == root:
              # This process has the requested element so extract it
              element[:] = mydata[myindex]

          # Broadcast from owner to everyone else
          comm.broadcast(element, root)

          # .. which is equivalent to ..

          if comm.rank == root:
              # We are root so send it to the other ranks
              for rank in range(comm.size):
                  if rank != root:
                      comm.send(element, rank, tag=123)
          else:
              # We don't have it so receive from root
              comm.receive(element, root, tag=123)

        """
        assert 0 <= root < self.size
        assert is_contiguous(a)
        self.comm.broadcast(a, root)

    def sendreceive(self, a, dest, b, src, sendtag=123, recvtag=123):
        assert 0 <= dest < self.size
        assert dest != self.rank
        assert is_contiguous(a)
        assert 0 <= src < self.size
        assert src != self.rank
        assert is_contiguous(b)
        return self.comm.sendreceive(a, dest, b, src, sendtag, recvtag)

    def send(self, a, dest, tag=123, block=True):
        assert 0 <= dest < self.size
        assert dest != self.rank
        assert is_contiguous(a)
        if not block:
            pass  # assert sys.getrefcount(a) > 3
        return self.comm.send(a, dest, tag, block)

    def ssend(self, a, dest, tag=123):
        assert 0 <= dest < self.size
        assert dest != self.rank
        assert is_contiguous(a)
        return self.comm.ssend(a, dest, tag)

    def receive(self, a, src, tag=123, block=True):
        assert 0 <= src < self.size
        assert src != self.rank
        assert is_contiguous(a)
        return self.comm.receive(a, src, tag, block)

    def test(self, request):
        """Test whether a non-blocking MPI operation has completed. A boolean
        is returned immediately and the request is not modified in any way.

        Parameters:

        request: MPI request
            Request e.g. returned from send/receive when block=False is used.

        """
        return self.comm.test(request)

    def testall(self, requests):
        """Test whether non-blocking MPI operations have completed. A boolean
        is returned immediately but requests may have been deallocated as a
        result, provided they have completed before or during this invokation.

        Parameters:

        request: MPI request
            Request e.g. returned from send/receive when block=False is used.

        """
        return self.comm.testall(requests)  # may deallocate requests!

    def wait(self, request):
        """Wait for a non-blocking MPI operation to complete before returning.

        Parameters:

        request: MPI request
            Request e.g. returned from send/receive when block=False is used.

        """
        self.comm.wait(request)

    def waitall(self, requests):
        """Wait for non-blocking MPI operations to complete before returning.

        Parameters:

        requests: list
            List of MPI requests e.g. aggregated from returned requests of
            multiple send/receive calls where block=False was used.

        """
        self.comm.waitall(requests)

    def abort(self, errcode):
        """Terminate MPI execution environment of all tasks in the group.
        This function only returns in the advent of an error occurring.

        Parameters:

        errcode: int
            Error code to return to the invoking environment.

        """
        return self.comm.abort(errcode)

    def name(self):
        """Return the name of the processor as a string."""
        return self.comm.name()

    def barrier(self):
        """Block execution until all process have reached this point."""
        self.comm.barrier()

    def get_members(self):
        """Return the subset of processes which are members of this MPI group
        in terms of the ranks they are assigned on the parent communicator.
        For the world communicator, this is all integers up to ``size``.

        Example::

          >>> world.rank, world.size
          (3, 4)
          >>> world.get_members()
          array([0, 1, 2, 3])
          >>> comm = world.new_communicator(array([2, 3]))
          >>> comm.rank, comm.size
          (1, 2)
          >>> comm.get_members()
          array([2, 3])
          >>> comm.get_members()[comm.rank] == world.rank
          True

        """
        return self.comm.get_members()

    def get_c_object(self):
        """Return the C-object wrapped by this debug interface.

        Whenever a communicator object is passed to C code, that object
        must be a proper C-object - *not* e.g. this debug wrapper.  For
        this reason.  The C-communicator object has a get_c_object()
        implementation which returns itself; thus, always call
        comm.get_c_object() and pass the resulting object to the C code.
        """
        c_obj = self.comm.get_c_object()
        assert type(c_obj) is _gpaw.Communicator
        return c_obj


# Serial communicator
class SerialCommunicator:
    size = 1
    rank = 0

    def __init__(self, parent=None):
        self.parent = parent

    def sum(self, array, root=-1):
        if isinstance(array, (int, float, complex)):
            return array

    def scatter(self, s, r, root):
        r[:] = s

    def min(self, value, root=-1):
        return value

    def max(self, value, root=-1):
        return value

    def broadcast(self, buf, root):
        pass

    def send(self, buff, root, tag=123, block=True):
        pass

    def barrier(self):
        pass

    def gather(self, a, root, b):
        b[:] = a

    def all_gather(self, a, b):
        b[:] = a

    def new_communicator(self, ranks):
        if self.rank not in ranks:
            return None
        return SerialCommunicator(parent=self)

    def test(self, request):
        return 1

    def testall(self, requests):
        return 1

    def wait(self, request):
        raise NotImplementedError('Calls to mpi wait should not happen in '
                                  'serial mode')

    def waitall(self, requests):
        if not requests:
            return
        raise NotImplementedError('Calls to mpi waitall should not happen in '
                                  'serial mode')

    def get_members(self):
        return np.array([0])

    def get_c_object(self):
        raise NotImplementedError('Should not get C-object for serial comm')


serial_comm = SerialCommunicator()

try:
    world = _gpaw.Communicator()
except AttributeError:
    world = serial_comm

    
class DryRunCommunicator(SerialCommunicator):
    def __init__(self, size=1, parent=None):
        self.size = size
        self.parent = parent
    
    def new_communicator(self, ranks):
        return DryRunCommunicator(len(ranks), parent=self)

    def get_c_object(self):
        return None  # won't actually be passed to C

if dry_run_size > 1:
    world = DryRunCommunicator(dry_run_size)


if debug:
    serial_comm = _Communicator(serial_comm)
    world = _Communicator(world)


size = world.size
rank = world.rank
parallel = (size > 1)


# XXXXXXXXXX for easier transition to Parallelization class
def distribute_cpus(parsize_domain, parsize_bands,
                    nspins, nibzkpts, comm=world,
                    idiotproof=True, mode='fd'):
    nsk = nspins * nibzkpts
    if mode in ['fd', 'lcao']:
        if parsize_bands is None:
            parsize_bands = 1
    else:
        # Plane wave mode:
        if parsize_bands is None:
            parsize_bands = comm.size // gcd(nsk, comm.size)

    p = Parallelization(comm, nsk)
    return p.build_communicators(domain=np.prod(parsize_domain),
                                 band=parsize_bands)

    
def old_distribute_cpus(parsize_domain, parsize_bands,
                        nspins, nibzkpts, comm=world,
                        idiotproof=True, mode='fd'):
    """Distribute k-points/spins to processors.

    Construct communicators for parallelization over
    k-points/spins and for parallelization using domain
    decomposition."""

    size = comm.size
    rank = comm.rank

    nsk = nspins * nibzkpts

    if mode in ['fd', 'lcao']:
        if parsize_bands is None:
            parsize_bands = 1

        if parsize_domain is not None:
            if type(parsize_domain) is int:
                ndomains = parsize_domain
            else:
                ndomains = (parsize_domain[0] *
                            parsize_domain[1] *
                            parsize_domain[2])
            assert (size // parsize_bands) % ndomains == 0

        else:
            ntot = nsk * parsize_bands
            ndomains = size // gcd(ntot, size)
    else:
        # Plane wave mode:
        ndomains = 1
        if parsize_bands is None:
            parsize_bands = size // gcd(nsk, size)

    assert size % parsize_bands == 0
        
    # How many spin/k-point combinations do we get per node:
    nu, x = divmod(nsk, size // parsize_bands // ndomains)
    assert x == 0 or nu >= 2 or not idiotproof, 'load imbalance!'

    r0 = (rank // ndomains) * ndomains
    ranks = np.arange(r0, r0 + ndomains)
    domain_comm = comm.new_communicator(ranks)

    r0 = rank % (ndomains * parsize_bands)
    ranks = np.arange(r0, r0 + size, ndomains * parsize_bands)
    kpt_comm = comm.new_communicator(ranks)

    r0 = rank % ndomains + kpt_comm.rank * (ndomains * parsize_bands)
    ranks = np.arange(r0, r0 + (ndomains * parsize_bands), ndomains)
    band_comm = comm.new_communicator(ranks)

    assert size == domain_comm.size * kpt_comm.size * band_comm.size

    return domain_comm, kpt_comm, band_comm


def compare_atoms(atoms, comm=world):
    """Check whether atoms objects are identical on all processors."""
    # Construct fingerprint:
    # ASE may return slightly different atomic positions (e.g. due
    # to MKL) so compare only first 8 decimals of positions
    fingerprint = np.array([md5_array(array, numeric=True) for array in
                            [atoms.positions.round(8),
                             atoms.cell,
                             atoms.pbc * 1.0,
                             atoms.get_initial_magnetic_moments()]])
    # Compare fingerprints:
    fingerprints = np.empty((comm.size, 4), fingerprint.dtype)
    comm.all_gather(fingerprint, fingerprints)
    mismatches = fingerprints.ptp(0)

    if debug:
        dumpfile = 'compare_atoms'
        for i in np.argwhere(mismatches).ravel():
            itemname = ['positions', 'cell', 'pbc', 'magmoms'][i]
            itemfps = fingerprints[:, i]
            itemdata = [atoms.positions,
                        atoms.cell,
                        atoms.pbc * 1.0,
                        atoms.get_initial_magnetic_moments()][i]
            if comm.rank == 0:
                print('DEBUG: compare_atoms failed for %s' % itemname)
                itemfps.dump('%s_fps_%s.pickle' % (dumpfile, itemname))
            itemdata.dump('%s_r%04d_%s.pickle' % (dumpfile, comm.rank,
                                                  itemname))

    # Use only the atomic positions from rank 0
    comm.broadcast(atoms.positions, 0)
    return not mismatches.any()

    
def broadcast(obj, root=0, comm=world):
    """Broadcast a Python object across an MPI communicator and return it."""
    if comm.rank == root:
        assert obj is not None
        string = pickle.dumps(obj, pickle.HIGHEST_PROTOCOL)
    else:
        assert obj is None
        string = None
    string = broadcast_string(string, root, comm)
    if comm.rank == root:
        return obj
    else:
        return pickle.loads(string)

        
def broadcast_string(string=None, root=0, comm=world):
    """Broadcast a Python string across an MPI communicator and return it.
    NB: Strings are immutable objects in Python, so the input is unchanged."""
    if comm.rank == root:
        assert isinstance(string, str)
        n = np.array(len(string), int)
    else:
        assert string is None
        n = np.zeros(1, int)
    comm.broadcast(n, root)
    if comm.rank == root:
        string = np.fromstring(string, np.int8)
    else:
        string = np.zeros(n, np.int8)
    comm.broadcast(string, root)
    return string.tostring()

    
def send_string(string, rank, comm=world):
    comm.send(np.array(len(string)), rank)
    comm.send(np.fromstring(string, np.int8), rank)

    
def receive_string(rank, comm=world):
    n = np.array(0)
    comm.receive(n, rank)
    string = np.empty(n, np.int8)
    comm.receive(string, rank)
    return string.tostring()

    
def alltoallv_string(send_dict, comm=world):
    scounts = np.zeros(comm.size, dtype=np.int)
    sdispls = np.zeros(comm.size, dtype=np.int)
    stotal = 0
    for proc in range(comm.size):
        if proc in send_dict:
            data = np.fromstring(send_dict[proc], np.int8)
            scounts[proc] = data.size
            sdispls[proc] = stotal
            stotal += scounts[proc]

    rcounts = np.zeros(comm.size, dtype=np.int)
    comm.alltoallv(scounts, np.ones(comm.size, dtype=np.int),
                   np.arange(comm.size, dtype=np.int),
                   rcounts, np.ones(comm.size, dtype=np.int),
                   np.arange(comm.size, dtype=np.int))
    rdispls = np.zeros(comm.size, dtype=np.int)
    rtotal = 0
    for proc in range(comm.size):
        rdispls[proc] = rtotal
        rtotal += rcounts[proc]
        rtotal += rcounts[proc]

    sbuffer = np.zeros(stotal, dtype=np.int8)
    for proc in range(comm.size):
        sbuffer[sdispls[proc]:(sdispls[proc] + scounts[proc])] = (
            np.fromstring(send_dict[proc], np.int8))

    rbuffer = np.zeros(rtotal, dtype=np.int8)
    comm.alltoallv(sbuffer, scounts, sdispls, rbuffer, rcounts, rdispls)

    rdict = {}
    for proc in range(comm.size):
        rdict[proc] = rbuffer[rdispls[proc]:(rdispls[proc] + rcounts[proc])].tostring()

    return rdict

    
def ibarrier(timeout=None, root=0, tag=123, comm=world):
    """Non-blocking barrier returning a list of requests to wait for.
    An optional time-out may be given, turning the call into a blocking
    barrier with an upper time limit, beyond which an exception is raised."""
    requests = []
    byte = np.ones(1, dtype=np.int8)
    if comm.rank == root:
        # Everybody else:
        for rank in range(0, root) + range(root + 1, comm.size):
            rbuf, sbuf = np.empty_like(byte), byte.copy()
            requests.append(comm.send(sbuf, rank, tag=2 * tag + 0,
                                      block=False))
            requests.append(comm.receive(rbuf, rank, tag=2 * tag + 1,
                                         block=False))
    else:
        rbuf, sbuf = np.empty_like(byte), byte
        requests.append(comm.receive(rbuf, root, tag=2 * tag + 0, block=False))
        requests.append(comm.send(sbuf, root, tag=2 * tag + 1, block=False))

    if comm.size == 1 or timeout is None:
        return requests

    t0 = time.time()
    while not comm.testall(requests):  # automatic clean-up upon success
        if time.time() - t0 > timeout:
            raise RuntimeError('MPI barrier timeout.')
    return []

    
def run(iterators):
    """Run through list of iterators one step at a time."""
    if not isinstance(iterators, list):
        # It's a single iterator - empty it:
        for i in iterators:
            pass
        return

    if len(iterators) == 0:
        return

    while True:
        try:
            results = [iter.next() for iter in iterators]
        except StopIteration:
            return results

            
class Parallelization:
    def __init__(self, comm, nspinkpts):
        self.comm = comm
        self.size = comm.size
        self.nspinkpts = nspinkpts
        
        self.kpt = None
        self.domain = None
        self.band = None
        
        self.nclaimed = 1
        self.navail = comm.size

    def set(self, kpt=None, domain=None, band=None):
        if kpt is not None:
            self.kpt = kpt
        if domain is not None:
            self.domain = domain
        if band is not None:
            self.band = band
        
        nclaimed = 1
        for group, name in zip([self.kpt, self.domain, self.band],
                               ['k-point', 'domain', 'band']):
            if group is not None:
                if self.size % group != 0:
                    msg = ('Cannot paralllize as the '
                           'communicator size %d is not divisible by the '
                           'requested number %d of ranks for %s '
                           'parallelization' % (self.size, group, name))
                    raise ValueError(msg)
                nclaimed *= group
        navail = self.size // nclaimed
        
        assert self.size % nclaimed == 0
        assert self.size % navail == 0

        self.navail = navail
        self.nclaimed = nclaimed

    def get_communicator_sizes(self, kpt=None, domain=None, band=None):
        self.set(kpt=kpt, domain=domain, band=band)
        self.autofinalize()
        return self.kpt, self.domain, self.band

    def build_communicators(self, kpt=None, domain=None, band=None,
                            order='kbd'):
        """Construct communicators.

        Returns a communicator for k-points, domains, bands and
        k-points/bands.  The last one "unites" all ranks that are
        responsible for the same domain.

        The order must be a permutation of the characters 'kbd', each
        corresponding to each a parallelization mode.  The last
        character signifies the communicator that will be assigned
        contiguous ranks, i.e. order='kbd' will yield contiguous
        domain ranks, whereas order='kdb' will yield contiguous band
        ranks."""
        self.set(kpt=kpt, domain=domain, band=band)
        self.autofinalize()
        
        comm = self.comm
        rank = comm.rank
        communicators = {}
        parent_stride = self.size
        offset = 0

        groups = dict(k=self.kpt, b=self.band, d=self.domain)

        # Build communicators in hierachical manner
        # The ranks in the first group have largest separation while
        # the ranks in the last group are next to each other
        for name in order:
            group = groups[name]
            stride = parent_stride // group
            # First rank in this group
            r0 = rank % stride + offset
            # Last rank in this group
            r1 = r0 + stride * group
            ranks = np.arange(r0, r1, stride)
            communicators[name] = comm.new_communicator(ranks)
            parent_stride = stride
            # Offset for the next communicator
            offset += communicators[name].rank * stride

        # We want a communicator for kpts/bands, i.e. the complement of the
        # grid comm: a communicator uniting all cores with the same domain.
        c1, c2, c3 = [communicators[name] for name in order]
        allranks = [range(c1.size), range(c2.size), range(c3.size)]
        
        def get_communicator_complement(name):
            relevant_ranks = list(allranks)
            relevant_ranks[order.find(name)] = [communicators[name].rank]
            ranks = np.array([r3 + c3.size * (r2 + c2.size * r1)
                              for r1 in relevant_ranks[0]
                              for r2 in relevant_ranks[1]
                              for r3 in relevant_ranks[2]])
            return comm.new_communicator(ranks)
        
        # The communicator of all processes that share a domain, i.e.
        # the combination of k-point and band dommunicators.
        communicators['D'] = get_communicator_complement('d')
        # For each k-point comm rank, a communicator of all
        # band/domain ranks.  This is typically used with ScaLAPACK
        # and LCAO orbital stuff.
        communicators['K'] = get_communicator_complement('k')
        return communicators
    
    def autofinalize(self):
        if self.kpt is None:
            self.set(kpt=self.get_optimal_kpt_parallelization())
        if self.domain is None:
            self.set(domain=self.navail)
        if self.band is None:
            self.set(band=self.navail)

        if self.navail > 1:
            raise RuntimeError('All the CPUs must be used')
    
    def get_optimal_kpt_parallelization(self, kptprioritypower=1.4):
        if self.domain and self.band:
            # Try to use all the CPUs for k-point parallelization
            ncpus = min(self.nspinkpts, self.navail)
            return ncpus
        ncpuvalues, wastevalues = self.find_kpt_parallelizations()
        scores = ((self.navail // ncpuvalues)
                  * ncpuvalues**kptprioritypower)**(1.0 - wastevalues)
        arg = np.argmax(scores)
        ncpus = ncpuvalues[arg]
        return ncpus

    def find_kpt_parallelizations(self):
        nspinkpts = self.nspinkpts
        ncpuvalues = []
        wastevalues = []
        
        ncpus = nspinkpts
        while ncpus > 0:
            if self.navail % ncpus == 0:
                nkptsmin = nspinkpts // ncpus
                nkptsmax = -(-nspinkpts // ncpus)
                effort = nkptsmax * ncpus
                efficiency = nspinkpts / float(effort)
                waste = 1.0 - efficiency
                wastevalues.append(waste)
                ncpuvalues.append(ncpus)
            ncpus -= 1
        return np.array(ncpuvalues), np.array(wastevalues)


def cleanup():
    error = getattr(sys, 'last_type', None)
    if error is not None:  # else: Python script completed or raise SystemExit
        if parallel and not (dry_run_size > 1):
            sys.stdout.flush()
            sys.stderr.write(('GPAW CLEANUP (node %d): %s occurred.  '
                              'Calling MPI_Abort!\n') % (world.rank, error))
            sys.stderr.flush()
            # Give other nodes a moment to crash by themselves (perhaps
            # producing helpful error messages)
            time.sleep(10)
            world.abort(42)

            
def exit(error='Manual exit'):
    # Note that exit must be called on *all* MPI tasks
    atexit._exithandlers = []  # not needed because we are intentially exiting
    if parallel and not (dry_run_size > 1):
        sys.stdout.flush()
        sys.stderr.write(('GPAW CLEANUP (node %d): %s occurred.  ' +
                          'Calling MPI_Finalize!\n') % (world.rank, error))
        sys.stderr.flush()
    else:
        cleanup(error)
    world.barrier()  # sync up before exiting
    sys.exit()  # quit for serial case, return to _gpaw.c for parallel case

atexit.register(cleanup)
