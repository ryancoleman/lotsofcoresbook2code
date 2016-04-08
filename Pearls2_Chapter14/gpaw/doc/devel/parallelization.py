#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from gpaw.mpi import world
from gpaw.utilities.dscftools import mpi_debug


W = world.size
N = 32

assert N%W == 0
M = N//W

# Create my share of data
data = np.arange(world.rank*M, (world.rank+1)*M)

# Let's calculate the global sum
slocal = data.sum()
s = world.sum(slocal)
mpi_debug('data: %s, slocal=%d, s=%d' % (data,slocal,s))
assert s == N*(N-1)//2

# Subtract the global mean
data -= s/N
mpi_debug('data: %s' % data)

# -------------------------------------------------------------------

if world.rank == 0:
    print('-'*16)

# Who has global index 11? The master needs it!
i = 11
rank, ilocal = divmod(i, M)
mpi_debug('rank=%d, ilocal=%d, i=%d' % (rank,ilocal,i))
assert rank*M + ilocal == i

# Do I have it?
if world.rank == rank:
    # Yes, so extract data (must be an array)
    idata = np.array([data[ilocal]], dtype=data.dtype)
else:
    # No, so just allocate space
    idata = np.empty(1, dtype=data.dtype)

# Broadcast from owner to everyone else
world.broadcast(idata, rank)

"""
# This does the same as broadcast with send/receive...

# Do I have it?
if world.rank == rank:
    # Yes, now send it to the others
    for other_rank in range(world.size):
        # We don't have to send it to ourselves
        if other_rank != rank:
            world.send(idata, other_rank, tag=123)
else:
    # No, so receive from the one that own the data
    world.receive(idata, rank, tag=123)
"""

mpi_debug('idata=%d' % idata)

# -------------------------------------------------------------------

if world.rank == 0:
    print('-'*16)

# The master just calculated auxilary data. Distribute it.
aux = np.empty(N, dtype=float)

# Only master knows the data right now
if world.rank == 0:
    np.random.seed(1234567)
    aux[:] = np.random.uniform(0,1,size=N).round(2)
    print('MASTER aux: %s, mean=%f' % (aux, aux.mean()))

# Allocate space for my part of the auxilary data
myaux = np.empty(M, dtype=float)

# Scatter parts from master to everyone
world.scatter(aux, myaux, 0)

"""
# This does the same as scatter with send/receive...

# Are we the master?
if world.rank == 0:
    # Yes, so extract my part directly
    myaux[:] = aux[0:M]

    # Now send parts to the slaves
    for slave_rank in range(1, world.size):
        youraux = aux[slave_rank*M:(slave_rank+1)*M]
        world.send(youraux, slave_rank, tag=123)
else:
    # No, so receive from the master
    world.receive(myaux, 0, tag=123)
"""

# We don't need original data anymore
del aux

# Try to calculate mean now
meanaux = world.sum(myaux.mean())/world.size
mpi_debug('myaux: %s, mean=%f' % (myaux,meanaux))

# -------------------------------------------------------------------

if world.rank == 0:
    print('-'*16)

# We've done something to our part of the auxilary data. Master needs it all
if world.rank == 0:
    result = np.empty(N, dtype=float)
else:
    result = None

# Do something to our auxilary data
myaux[:] = np.sin(2*np.pi*myaux).round(3)
mpi_debug('myaux: %s' % myaux)

# Gather parts from everyone on the master
world.gather(myaux, 0, result)

"""
# This does the same as gather with send/receive...

# Are we the master?
if world.rank == 0:
    # Yes, so extract my part directly
    result[0:M] = myaux[:]

    # Now receive parts from the slaves
    for slave_rank in range(1, world.size):
        youraux = np.empty(M, dtype=float)
        world.receive(youraux, slave_rank, tag=123)
        result[slave_rank*M:(slave_rank+1)*M] = youraux
else:
    # No, so send to the master
    world.send(myaux, 0, tag=123)
"""

mpi_debug('result: %s' % result)

