from __future__ import print_function
from gpaw.mpi import world
from gpaw.blacs import BlacsGrid, Redistributor

if world.size < 2:
    raise ValueError('Runs on two or more processors')

grid = BlacsGrid(world, 2, world.size // 2)

desc = grid.new_descriptor(12, 8, 2, 3)

a = desc.zeros()
a[:] = world.rank

subdesc = grid.new_descriptor(7, 7, 2, 2)
b = subdesc.zeros()

r = Redistributor(grid.comm, desc, subdesc)

ia = 3
ja = 2
ib = 1
jb = 1
M = 4
N = 5

r.redistribute(a, b, M, N, ia, ja, ib, jb)

a0 = desc.collect_on_master(a)
b0 = subdesc.collect_on_master(b)
if world.rank == 0:
    print(a0)
    print(b0)
    xa = a0[ia:ia + M, ja:ja + N]
    xb = b0[ib:ib + M, jb:jb + N]
    assert (xa == xb).all()
    
