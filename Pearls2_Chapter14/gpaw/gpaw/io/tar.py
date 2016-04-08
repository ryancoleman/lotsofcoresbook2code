import os
import time
import tarfile
import xml.sax

import numpy as np

from gpaw.mpi import broadcast as mpi_broadcast
from gpaw.mpi import world

from gpaw.io import FileReference


intsize = 4
floatsize = np.array([1], float).itemsize
complexsize = np.array([1], complex).itemsize
itemsizes = {'int': intsize, 'float': floatsize, 'complex': complexsize}

    
class Writer:
    def __init__(self, name, comm=world):
        self.comm = comm # for possible future use
        self.dims = {}
        self.files = {}
        self.xml1 = ['<gpaw_io version="0.1" endianness="%s">' %
                     ('big', 'little')[int(np.little_endian)]]
        self.xml2 = []
        if os.path.isfile(name):
            os.rename(name, name[:-4] + '.old'+name[-4:])
        self.tar = tarfile.open(name, 'w')
        self.mtime = int(time.time())
        
    def dimension(self, name, value):
        if name in self.dims.keys() and self.dims[name] != value:
            raise Warning('Dimension %s changed from %s to %s' % \
                          (name, self.dims[name], value))
        self.dims[name] = value

    def __setitem__(self, name, value):
        if isinstance(value, float):
            value = repr(value)
        self.xml1 += ['  <parameter %-20s value="%s"/>' %
                      ('name="%s"' % name, value)]
        
    def add(self, name, shape, array=None, dtype=None, units=None,
            parallel=False, write=True):
        if array is not None:
            array = np.asarray(array)

        self.dtype, type, itemsize = self.get_data_type(array, dtype)
        self.xml2 += ['  <array name="%s" type="%s">' % (name, type)]
        self.xml2 += ['    <dimension length="%s" name="%s"/>' %
                      (self.dims[dim], dim)
                      for dim in shape]
        self.xml2 += ['  </array>']
        self.shape = [self.dims[dim] for dim in shape]
        size = itemsize * int(np.product([self.dims[dim] for dim in shape]))
        self.write_header(name, size)
        if array is not None:
            self.fill(array)

    def get_data_type(self, array=None, dtype=None):
        if dtype is None:
            dtype = array.dtype

        if dtype in [int, bool]:
            dtype = np.int32

        dtype = np.dtype(dtype)
        type = {np.int32: 'int',
                np.float64: 'float',
                np.complex128: 'complex'}[dtype.type]

        return dtype, type, dtype.itemsize

    def fill(self, array, *indices, **kwargs):
        self.write(np.asarray(array, self.dtype).tostring())

    def write_header(self, name, size):
        assert name not in self.files.keys()
        tarinfo = tarfile.TarInfo(name)
        tarinfo.mtime = self.mtime
        tarinfo.size = size
        self.files[name] = tarinfo
        self.size = size
        self.n = 0
        self.tar.addfile(tarinfo)

    def write(self, string):
        self.tar.fileobj.write(string)
        self.n += len(string)
        if self.n == self.size:
            blocks, remainder = divmod(self.size, tarfile.BLOCKSIZE)
            if remainder > 0:
                self.tar.fileobj.write('\0' * (tarfile.BLOCKSIZE - remainder))
                blocks += 1
            self.tar.offset += blocks * tarfile.BLOCKSIZE
        
    def close(self):
        self.xml2 += ['</gpaw_io>\n']
        string = '\n'.join(self.xml1 + self.xml2)
        self.write_header('info.xml', len(string))
        self.write(string)
        self.tar.close()


class Reader(xml.sax.handler.ContentHandler):
    def __init__(self, name, comm=world):
        self.comm = comm # used for broadcasting replicated data
        self.master = (self.comm.rank == 0)
        self.dims = {}
        self.shapes = {}
        self.dtypes = {}
        self.parameters = {}
        xml.sax.handler.ContentHandler.__init__(self)
        self.tar = tarfile.open(name, 'r')
        f = self.tar.extractfile('info.xml')
        xml.sax.parse(f, self)

    def startElement(self, tag, attrs):
        if tag == 'gpaw_io':
            self.byteswap = ((attrs['endianness'] == 'little')
                             != np.little_endian)
        elif tag == 'array':
            name = attrs['name']
            self.dtypes[name] = attrs['type']
            self.shapes[name] = []
            self.name = name
        elif tag == 'dimension':
            n = int(attrs['length'])
            self.shapes[self.name].append(n)
            self.dims[attrs['name']] = n
        else:
            assert tag == 'parameter'
            try:
                value = eval(attrs['value'], {})
            except (SyntaxError, NameError):
                value = attrs['value'].encode()
            self.parameters[attrs['name']] = value

    def dimension(self, name):
        return self.dims[name]
    
    def __getitem__(self, name):
        return self.parameters[name]

    def has_array(self, name):
        return name in self.shapes
    
    def get(self, name, *indices, **kwargs):
        broadcast = kwargs.pop('broadcast', False)
        if self.master or not broadcast:
            fileobj, shape, size, dtype = self.get_file_object(name, indices)
            array = np.fromstring(fileobj.read(size), dtype)
            if self.byteswap:
                array = array.byteswap()
            if dtype == np.int32:
                array = np.asarray(array, int)
            array.shape = shape
            if shape == ():
                array = array.item()
        else:
            array = None

        if broadcast:
            array = mpi_broadcast(array, 0, self.comm)
        return array
    
    def get_reference(self, name, indices, length=None):
        fileobj, shape, size, dtype = self.get_file_object(name, indices)
        assert dtype != np.int32
        return TarFileReference(fileobj, shape, dtype, self.byteswap, length)
    
    def get_file_object(self, name, indices):
        dtype, type, itemsize = self.get_data_type(name)
        fileobj = self.tar.extractfile(name)
        n = len(indices)
        shape = self.shapes[name]
        size = itemsize * np.prod(shape[n:], dtype=int)
        offset = 0
        stride = size
        for i in range(n - 1, -1, -1):
            offset += indices[i] * stride
            stride *= shape[i]
        fileobj.seek(offset)
        return fileobj, shape[n:], size, dtype

    def get_data_type(self, name):
        type = self.dtypes[name]
        dtype = np.dtype({'int': np.int32,
                          'float': float,
                          'complex': complex}[type])
        return dtype, type, dtype.itemsize

    def get_parameters(self):
        return self.parameters

    def close(self):
        self.tar.close()


class TarFileReference(FileReference):
    def __init__(self, fileobj, shape, dtype, byteswap, length):
        self.fileobj = fileobj
        self.shape = tuple(shape)
        self.dtype = dtype
        self.itemsize = dtype.itemsize
        self.byteswap = byteswap
        self.offset = fileobj.tell()
        self.length = length

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, indices):
        if isinstance(indices, slice):
            start, stop, step = indices.indices(len(self))
            if start != 0 or step != 1 or stop != len(self):
                raise NotImplementedError('You can only slice a TarReference '
                                          'with [:] or [int]')
            else:
                indices = ()
        elif isinstance(indices, int):
            indices = (indices,)
        else: # Probably tuple or ellipsis
            raise NotImplementedError('You can only slice a TarReference '
                                      'with [:] or [int]')
            
        n = len(indices)

        size = np.prod(self.shape[n:], dtype=int) * self.itemsize
        offset = self.offset
        stride = size
        for i in range(n - 1, -1, -1):
            offset += indices[i] * stride
            stride *= self.shape[i]
        self.fileobj.seek(offset)
        array = np.fromstring(self.fileobj.read(size), self.dtype)
        if self.byteswap:
            array = array.byteswap()
        array.shape = self.shape[n:]
        if self.length:
            array = array[..., :self.length].copy()
        return array
