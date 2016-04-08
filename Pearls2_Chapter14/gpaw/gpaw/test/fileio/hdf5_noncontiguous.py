import numpy as np
import os
from gpaw.io.hdf5_highlevel import File
from gpaw.mpi import world, rank

# write a slice of NumPy array
data = np.arange(10, dtype=float)
sub = data[::2]

if world.size > 1:
    comm = world.get_c_object()
else:
    comm = None

file = File('tmp.hdf5', 'w', comm=comm)

dset = file.create_dataset('noncont', sub.shape, sub.dtype) 
if rank == 0:
    selection = 'all'
else:
    selection = None
dset.write(sub, selection)
dset.close()
file.close()

# read data back
file = File('tmp.hdf5', 'r', comm=comm)
dset = file['noncont']
new_data = np.ndarray(dset.shape, dset.dtype, order='C')
dset.read(new_data)
assert((new_data == sub).all)
# read also to noncontiguous array
new_data = np.zeros(10, float)
new_sub = new_data[::2]
dset.read(new_sub)
assert((new_sub == sub).all)

dset.close()
file.close()
if rank == 0:
    os.remove('tmp.hdf5')
