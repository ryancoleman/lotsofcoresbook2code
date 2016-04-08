"""Light weight Python interface to HDF5""" 

#  Copyright (C) 2010-2011     CSC - IT Center for Science Ltd.
#  Please see the accompanying LICENSE file for further information.

# The module is heavily inspired by h5py, however, no 
# h5py code is used directly except where explicitly noted

import numpy as np
from _gpaw_hdf5 import *

def numpy_type_from_h5(datatype):
    """Simple conversion from HDF5 datatype to NumPy dtype"""
    cls = h5t_get_class(datatype)
    if cls == H5T_INTEGER:
        dtype = int
    elif cls == H5T_FLOAT:
        dtype = float
    elif cls == H5T_COMPOUND:
        dtype = complex
    elif cls == H5T_ENUM:
        dtype = bool
    elif cls == H5T_STRING:
        str_size = h5t_get_size(datatype)
        dtype = np.dtype('|S' + str(str_size))
    else:
        raise RuntimeError('Unsupported HDF5 datatype')

    return dtype

class Iterator:
    """Iterator over the datasets and subgroups in a File or in a Group."""

    def __init__(self, obj):
        self.idx = 0
        self.obj = obj
        self.nobj = obj.get_num_objs()

    def __iter__(self):
        return self

    def next(self):
        if self.idx == self.nobj:
            self.obj = None
            raise StopIteration

        key = self.obj.get_name_by_idx(self.idx)
        self.idx += 1
        return key


class Group:
    """This class defines a HDF5 group."""
    def __init__(self, loc_id, name, create=False):
        if create:
            self.id = h5g_create(loc_id, name)
        else:
            self.id = h5g_open(loc_id, name)

        self.attrs = Attributes(self.id)

        self.name = name
        self.opened = True

    def create_group(self, name):
        return Group(self.id, name, create=True)

    def create_dataset(self, name, shape, dtype):
        """Create a dataset with the NumPy equivalent type  and shape.

           Arguments:
           ----------
           dtype : Python or NumPy type
           shape : tuple containing the shape of dataset
           """

        return Dataset(self.id, name, dtype, shape, create=True)

    def get_num_objs(self):
        return h5g_get_num_objs(self.id)

    def get_name_by_idx(self, idx):
        return h5l_get_name_by_idx(self.id, idx)

    def close(self):
        h5g_close(self.id)
        self.opened = False

    # Dictionary like interface
    def keys(self):
        return list(self)

    def values(self):
        return [self[x] for x in self.keys()]

    def items(self):
        return [(x, self[x]) for x in self.keys()]

    def __getitem__(self, key):
        oid = h5o_open(self.id, key)
        obj_type = h5i_get_type(oid)
        h5o_close(oid)
        if obj_type == H5I_GROUP:
            return Group(self.id, key)
        elif obj_type == H5I_DATASET:
            return Dataset(self.id, key)
        else:
           # This shoul never happen
           raise RuntimeError('Accessing unknown object type')

    def __iter__(self):
        return Iterator(self)

    def __len__(self):
        return self.get_num_objs()

    def __del__(self):
        if self.opened:
            self.close()

class File(Group):
    """This class defines a HDF5 file."""

    def __init__(self, name, mode='r', comm=None):
        """If a communicator is passed as an argument, file is 
           created or opened for parallel IO."""

        if comm is None:
            plist = H5P_DEFAULT
        else:
            plist = h5p_create(H5P_FILE_ACCESS)
            h5p_set_fapl_mpio(plist, comm)

        if mode == 'r':
            self.id = h5f_open(name, mode, plist)
        elif mode == 'w':
            self.id = h5f_create(name, plist)
        else:
            raise RuntimeError('Unsupported file open/create mode')

        self.attrs = Attributes(self.id)

        if comm is not None:
            h5p_close(plist)

        self.name = name
        self.opened = True

    def close(self):
        h5f_close(self.id)
        self.opened = False

class Dataset(object):
    """This class defines a HDF5 dataset.

       Attributes:

    ============= =============================================
    Name           Description
    ============= =============================================
    ``shape``     Tuple describing the dataspace of dataset
    ``dtype``     NumPy equivalent datatype of dataset
    ``dataspace`` HDF5 dataspace identifier (integer)
    ``datatype``  HDF5 datatype identifier (integer)
    ``id``        HDF5 identifier of dataset (integer)
    ============= =============================================
"""

    def __init__(self, loc_id, name, dtype=None, shape=None, create=False):
        if create:
            self.shape = shape
            self.dtype = dtype
            self.dataspace = h5s_create(np.asarray(shape))
            self.datatype = h5_type_from_numpy(np.ndarray((1,), dtype))
            self.id = h5d_create(loc_id, name, self.datatype, self.dataspace)
        else:
            self.id = h5d_open(loc_id, name)
            self.dataspace = h5d_get_space(self.id)
            self.datatype = h5d_get_type(self.id)
            self.shape = h5s_get_shape(self.dataspace)
            self.dtype = numpy_type_from_h5(self.datatype)

        self.attrs = Attributes(self.id)

        self.name = name
        self.opened = True

    def write(self, data, selection='all', collective=False):
        """Write NumPy array data into the dataset.

           selection can be a hyperslab instance for writing to a part
           of the dataset or None for no reading 
           collective specifies the use of collective IO."""

        if collective:
            plist = h5p_create(H5P_DATASET_XFER)
            h5p_set_dxpl_mpio(plist)
        else:
            plist = H5P_DEFAULT

        filespace = self.dataspace
        memspace = h5s_create(np.asarray(data.shape))
        memtype = h5_type_from_numpy(np.ndarray((1,), data.dtype))

        if selection is None:
            h5s_select_none(memspace)
            h5s_select_none(filespace)

        if isinstance(selection, HyperslabSelection):
            selection.select(self)

        # data array has to be contiguous, if not make a copy
        if not data.flags.contiguous:
            data = data.copy()
        h5d_write(self.id, memtype, memspace, filespace, data, plist)
        
        h5s_close(memspace)
        h5t_close(memtype)
        if collective:
            h5p_close(plist)

    def read(self, data, selection='all', collective=False):
        """Read the dataset into NumPy array data

           selection can be a hyperslab instance for reading part 
           of the dataset or None for no reading 
           collective specifies the use of collective IO."""

        if collective:
            plist = h5p_create(H5P_DATASET_XFER)
            h5p_set_dxpl_mpio(plist)
        else:
            plist = H5P_DEFAULT

        filespace = self.dataspace
        memspace = h5s_create(np.asarray(data.shape))
        memtype = h5_type_from_numpy(np.ndarray((1,), data.dtype))

        if selection is None:
            h5s_select_none(memspace)
            h5s_select_none(filespace)

        if isinstance(selection, HyperslabSelection):
            selection.select(self)

        # data array has to be contiguous, if not, create a 
        # temporary array 
        if not data.flags.contiguous:
            data_buf = data.copy()
        else:
            data_buf = data 
        h5d_read(self.id, memtype, memspace, filespace, data_buf, plist)

        if not data.flags.contiguous:
            data[:] = data_buf[:]
        
        h5s_close(memspace)
        h5t_close(memtype)
        if collective:
            h5p_close(plist)

    # Interface for access with slicing syntax
    def __setitem__(self, indices, value):
        sel = HyperslabSelection(indices, self.shape)

        if not isinstance(value, np.ndarray):
            value = np.asarray(value)
        if sel.mshape != value.shape:
            raise RuntimeError('Shapes not compatible')
        self.write(value, selection=sel)

    def __getitem__(self, indices):
        sel = HyperslabSelection(indices, self.shape)
        array = np.zeros(sel.mshape, self.dtype)
        self.read(array, sel)
        return array

    def close(self):
        h5t_close(self.datatype)
        h5s_close(self.dataspace)
        h5d_close(self.id)
        self.opened = False
        
    def __del__(self):
        if self.opened:
            self.close()

class Attributes:
    """A dictionary like interface to HDF5 attributes.

       Attributes can be written with the::

          attrs['name'] = value

       and read with the::

          value = attrs['name'] 

       syntax. 
       Values are returned always as NumPy arrays.

    """

    def __init__(self, parent_id):
        self.loc_id = parent_id

    def __setitem__(self, key, data):
        # Should we delete existing attributes?
        # Special treatment for possible None data
        if data is None:
            data = repr(data)

        data = np.asarray(data)
        dataspace = h5s_create(np.asarray(data.shape))
        datatype = h5_type_from_numpy(np.ndarray((1, ), data.dtype))
        id = h5a_create(self.loc_id, key, datatype, dataspace)
        # ensure that data is contiguous
        data = data.copy()
        h5a_write(id, datatype, data)
        h5s_close(dataspace)
        h5t_close(datatype)
        h5a_close(id)

    def __getitem__(self, key):
        id = h5a_open(self.loc_id, key)
        dataspace = h5a_get_space(id)
        datatype = h5a_get_type(id)
        shape = h5s_get_shape(dataspace)
        dtype = numpy_type_from_h5(datatype)
        memtype = h5_type_from_numpy(np.ndarray((1,), dtype))
        data = np.empty(shape, dtype)
        h5a_read(id, memtype, data)
        h5s_close(dataspace)
        h5t_close(datatype)
        h5t_close(memtype)
        h5a_close(id)

        if len(shape) == 0:
            data = np.asscalar(data)
        return data

    def __contains__(self, name):
        return h5a_exists_by_name(self.loc_id, name)

# The following four functions 
# _expand_ellipsis, _translate_int, _translate_slice and _handle_simple
# are direct copy-paste from h5py
# Copyright (C) 2008 Andrew Collette
# http://h5py.alfven.org
# License: New BSD
def _expand_ellipsis(args, rank):
    """ Expand ellipsis objects and fill in missing axes.
    """
    n_el = sum(1 for arg in args if arg is Ellipsis)
    if n_el > 1:
        raise ValueError("Only one ellipsis may be used.")
    elif n_el == 0 and len(args) != rank:
        args = args + (Ellipsis,)

    final_args = []
    n_args = len(args)
    for idx, arg in enumerate(args):
        if arg is Ellipsis:
            final_args.extend( (slice(None,None,None),)*(rank-n_args+1) )
        else:
            final_args.append(arg)

    if len(final_args) > rank:
        raise TypeError("Argument sequence too long")

    return final_args

def _translate_int(exp, length):
    """ Given an integer index, return a 3-tuple
        (start, count, step)
        for hyperslab selection
    """
    if exp < 0:
        exp = length+exp

    if not 0<=exp<length:
        raise ValueError("Index (%s) out of range (0-%s)" % (exp, length-1))
    return exp, 1, 1

def _translate_slice(exp, length):
    """ Given a slice object, return a 3-tuple
        (start, count, step)
        for use with the hyperslab selection routines
    """
    start, stop, step = exp.start, exp.stop, exp.step
    if start is None:
        start = 0
    else:
        start = int(start)
    if stop is None:
        stop = length
    else:
        stop = int(stop)
    if step is None:
        step = 1
    else:
        step = int(step)
    if step < 1:
        raise ValueError("Step must be >= 1 (got %d)" % step)
    if stop == start:
        raise ValueError("Zero-length selections are not allowed")
    if stop < start:
        raise ValueError("Reverse-order selections are not allowed")
    if start < 0:
        start = length+start
    if stop < 0:
        stop = length+stop

    if not 0 <= start <= (length-1):
        raise ValueError("Start index %s out of range (0-%d)" % (start,
length-1))
    if not 1 <= stop <= length:
        raise ValueError("Stop index %s out of range (1-%d)" % (stop, length))
    count = (stop-start)//step
    if (stop-start) % step != 0:
        count += 1
    
    if start+count > length:
        raise ValueError("Selection out of bounds (%d; axis has %d)" %
(start+count,length))
    
    return start, count, step
    
def _handle_simple(shape, args):
    """ Process a "simple" selection tuple, containing only slices and
        integer objects.  Return is a 4-tuple with tuples for start,
        count, step, and a flag which tells if the axis is a "scalar"
        selection (indexed by an integer).
    
        If "args" is shorter than "shape", the remaining axes are fully
        selected.
    """
    args = _expand_ellipsis(args, len(shape))

    start = []
    count = []
    step  = []
    scalar = []
    
    for arg, length in zip(args, shape):
        if isinstance(arg, slice):
            x,y,z = _translate_slice(arg, length)
            s = False
        else:
            try:
                x,y,z = _translate_int(int(arg), length)
                s = True
            except TypeError:
                raise TypeError('Illegal index "%s" (must be a slice or number)' % arg)
        start.append(x)
        count.append(y)
        step.append(z)
        scalar.append(s)

    return tuple(start), tuple(count), tuple(step), tuple(scalar)

class HyperslabSelection:
    """This class defines a simple hyperslab selection for a HDF5 dataspace.

    In HDF5, hyperslab is defined by offset, stride, count, and block 
    arguments, in the Python side simple indeces or slices are used
    (start, stop, step). The block argument of HDF5 is not used here. """
   
    def __init__(self, indices, shape):
        """Create a hyperslab descriptor."""
        if not isinstance(indices, tuple):
            indices = (indices,)

        dest_shape = shape
        self.indices = indices

        start, count, step, scalar = _handle_simple(dest_shape, indices)
        # This class follows HDF5 convention offset, stride, count
        self.hyperslab = (np.array(start), np.array(step), np.array(count))
        self.mshape = tuple(x for x, y in zip(count, scalar) if not y)

    def select(self, dataset):
        """Make the actual hyperslab selection to the dataspace of dataset."""
        dataspace = dataset.dataspace
        offset, stride, count = self.hyperslab
        h5s_select_hyperslab(dataspace, offset, stride, count, None)

