from __future__ import print_function
from gpaw.fd_operators import Gradient
import numpy as np
from gpaw.grid_descriptor import GridDescriptor
from gpaw.mpi import world

if world.size > 4:
    # Grid is so small that domain decomposition cannot exceed 4 domains
    assert world.size % 4 == 0
    group, other = divmod(world.rank, 4)
    ranks = np.arange(4*group, 4*(group+1))
    domain_comm = world.new_communicator(ranks)
else:
    domain_comm = world

gd = GridDescriptor((8, 1, 1), (8.0, 1.0, 1.0), comm=domain_comm)
a = gd.zeros()
dadx = gd.zeros()
a[:, 0, 0] = np.arange(gd.beg_c[0], gd.end_c[0])
gradx = Gradient(gd, v=0)
print(a.itemsize, a.dtype, a.shape)
print(dadx.itemsize, dadx.dtype, dadx.shape)
gradx.apply(a, dadx)

#   a = [ 0.  1.  2.  3.  4.  5.  6.  7.]
#
#   da
#   -- = [-2.5  1.   1.   1.   1.   1.  1.  -2.5]
#   dx

dadx = gd.collect(dadx, broadcast=True)
assert dadx[3, 0, 0] == 1.0 and np.sum(dadx[:, 0, 0]) == 0.0

gd = GridDescriptor((1, 8, 1), (1.0, 8.0, 1.0), (1, 0, 1), comm=domain_comm)
dady = gd.zeros()
a = gd.zeros()
grady = Gradient(gd, v=1)
a[0, :, 0] = np.arange(gd.beg_c[1], gd.end_c[1]) - 1
grady.apply(a, dady)

#   da
#   -- = [0.5  1.   1.   1.   1.   1.  -2.5]
#   dy

dady = gd.collect(dady, broadcast=True)
assert dady[0, 0, 0] == 0.5 and np.sum(dady[0, :, 0]) == 3.0

# a GUC cell
gd = GridDescriptor((1, 7, 1),
                    ((1.0, 0.0, 0.0),
                     (5.0, 5.0, 0.0),
                     (0.0, 0.0, 0.7)), comm=domain_comm)
dady = gd.zeros()
grady = Gradient(gd, v=1)
a = gd.zeros()
a[0, :, 0] = np.arange(gd.beg_c[1], gd.end_c[1]) - 1
grady.apply(a, dady)

#   da
#   -- = [-3.5  1.4   1.4   1.4   1.4   1.4  -3.5]
#   dy

dady = gd.collect(dady, broadcast=True)
assert dady[0, 0, 0] == -3.5 and abs(np.sum(dady[0, :, 0])) < 1E-12

