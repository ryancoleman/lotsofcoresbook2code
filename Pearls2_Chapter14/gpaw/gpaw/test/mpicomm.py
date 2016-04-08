
import numpy as np
from gpaw import debug, dry_run
from gpaw.mpi import world, serial_comm, _Communicator, \
                     SerialCommunicator, DryRunCommunicator

even_comm = world.new_communicator(np.arange(0, world.size, 2))
if world.size > 1:
    odd_comm = world.new_communicator(np.arange(1, world.size, 2))
else:
    odd_comm = None

if world.rank % 2 == 0:
    assert odd_comm is None
    comm = even_comm
else:
    assert even_comm is None
    comm = odd_comm

hasmpi = False
try:
    import _gpaw
    hasmpi = hasattr(_gpaw, 'Communicator')
except ImportError, AttributeError:
    pass

assert world.parent is None
assert comm.parent is world
if hasmpi:
    assert comm.parent.get_c_object() is world.get_c_object()
    assert comm.get_c_object().parent is world.get_c_object()

commranks = np.arange(world.rank % 2, world.size, 2)
assert np.all(comm.get_members() == commranks)
assert comm.get_members()[comm.rank] == world.rank

subcomm = comm.new_communicator(np.array([comm.rank]))
assert subcomm.parent is comm
assert subcomm.rank == 0 and subcomm.size == 1
assert subcomm.get_members().item() == comm.rank

if debug:
    assert isinstance(world, _Communicator)
    assert isinstance(comm, _Communicator)
    assert isinstance(subcomm, _Communicator)
elif world is serial_comm:
    assert isinstance(world, SerialCommunicator)
    assert isinstance(comm, SerialCommunicator)
    assert isinstance(subcomm, SerialCommunicator)
elif hasmpi:
    assert isinstance(world, _gpaw.Communicator)
    assert isinstance(comm, _gpaw.Communicator)
    assert isinstance(subcomm, _gpaw.Communicator)
