import numpy as np
import os
from gpaw.io.hdf5_highlevel import File, HyperslabSelection
from gpaw.mpi import world, rank, size

# data to be written
attrs = {'int' : 42, 'float' : 6.2, 'complex' : 1.0 + 2.0j, 
         'bool1' : False, 'bool2' : True, 'data' : np.array(((1,2),(3,4)))}
datasize = 4
data = {}
for dtype in ['int', 'float', 'complex']:
    data[dtype] = np.arange(rank * datasize, rank * datasize + datasize, 
                            dtype = dtype)
data['bool'] = np.array([True, False, False, True])

def write(attrs, data, file, parallel):
    file_title = 'test file'
    file.attrs['title'] = file_title
    grp = file.create_group('group1')
    # write attributes
    for key, val in attrs.items():
        grp.attrs[key] = val
    # write datasets
    shape = (datasize * size,)
    for key, val in data.items():
        dset = grp.create_dataset(key, shape, val.dtype) 
        indices = slice(rank * datasize, (rank + 1) * datasize)
        selection = HyperslabSelection(indices, dset.shape)
        dset.write(val, selection, collective=parallel)

def read(attrs, data, file, parallel):
    file_title = file.attrs['title']
    assert(file_title == 'test file')
    grp = file['group1']
    # read attributes
    for key, val in attrs.items():
        new_val = grp.attrs[key]
        if key == 'data':
            assert((new_val == val).all())
        else:
            assert(new_val == val)
    # read data
    for key, val in data.items():
        dset = grp[key]
        indices = slice((size - rank - 1) * datasize, 
                        (size - rank) * datasize)
        selection = HyperslabSelection(indices, dset.shape)
        new_val = np.ndarray(selection.mshape, dset.dtype, order='C')
        dset.read(new_val, selection, collective=parallel)
        assert((new_val == val).all)
    # read data in reversed order
    for key in data:
        if key == 'bool':
            continue
        dset = grp[key]
        indices = slice((size - rank - 1) * datasize, 
                        (size - rank) * datasize)
        selection = HyperslabSelection(indices, dset.shape)
        new_val = np.ndarray(selection.mshape, dset.dtype, order='C')
        dset.read(new_val, selection, collective=parallel)
        ref_val = np.arange((size - rank - 1) * datasize, 
                            (size - rank) * datasize, dtype=key)
        assert((new_val == ref_val).all)

if world.size > 1:
    comm = world.get_c_object()
    parallel = True
else:
    comm = None
    parallel = False

file = File('tmp.hdf5', 'w', comm=comm)
write(attrs, data, file, parallel)
file.close()

file = File('tmp.hdf5', 'r', comm=comm)
read(attrs, data, file, parallel)
file.close()

if rank == 0:
    os.remove('tmp.hdf5')
