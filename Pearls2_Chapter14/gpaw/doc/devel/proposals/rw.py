from __future__ import print_function
"""

File content::
    
    0: "IOASE..."
    8: version
    16: nitems (int64)
    24: 32 (position of offsets, int64)
    32: p0 (offset to json data, int64)
    40: 8-byte aligned ndarrays
    p0: n (length of json data, int64)
    p0+8: json data
    p0+8+n: EOF

"""

# magig prefix?, ascii header? See hdf5 header,
# ordereddict, endianness, todict?

import numpy as np

from ase.db.jsondb import encode, decode


VERSION = 1
N1 = 42  # block size - max number of items: 1, N1, N1*N1, N1*N1*N1, ...


def align(fd):
    """Advance file descriptor to 8 byte alignment and return position."""
    pos = fd.tell()
    r = pos % 8
    if r == 0:
        return pos
    fd.write(b'#' * (8 - r))
    return pos + 8 - r


def writeint(fd, n, pos=None):
    """Write 64 bit integer n at pos or current position."""
    if pos is not None:
        fd.seek(pos)
    np.array(n, np.int64).tofile(fd)
    

class Writer:
    def __init__(self, fd, mode='w', data=None):
        """Create writer object.

        The data dictionary holds:

        * data for type bool, int, float, complex and str
        * shape and dtype for ndarrays
        * class names for other objects

        These other objects must have a write() method and a static
        read() method."""

        assert np.little_endian
        
        if data is None:
            data = {}
            if mode == 'w':
                self.nitems = 0
                self.itemoffsets = 32
                self.offsets = np.array([-1], np.int64)

                fd = open(fd, 'wb')
            
                # Write file format identifier:
                fd.write(b'IOASE...')
                np.array([VERSION, self.nitems, self.itemoffsets],
                         np.int64).tofile(fd)
                self.offsets.tofile(fd)
            elif mode == 'a':
                fd = open(fd, 'r+b')
            
                version, self.nitems, self.itemoffsets, offsets = \
                    read_header(fd)
                assert version == VERSION
                n = 1
                while self.nitems > n:
                    n *= N1
                padding = np.zeros(n - self.nitems, np.int64)
                self.offsets = np.concatenate((offsets, padding))
                fd.seek(0, 2)
            else:
                2 / 0
            
        self.fd = fd
        self.data = data
        
        # Shape and dtype of array beeing filled:
        self.shape = (0,)
        self.dtype = None
        
    def add_array(self, name, shape, dtype=float, delayed_read=True):
        if isinstance(shape, int):
            shape = (shape,)
            
        i = align(self.fd)
        self.data[name] = {'_type': 'numpy.ndarray',
                           'shape': shape,
                           'dtype': np.dtype(dtype).name,
                           'offset': i}
            
        if delayed_read:
            self.data[name]['_delayed'] = True
            
        assert self.shape[0] == 0, 'last array not done'
        
        self.dtype = dtype
        self.shape = shape
        
    def fill(self, a):
        assert a.dtype == self.dtype
        if a.shape[1:] == self.shape[1:]:
            assert a.shape[0] <= self.shape[0]
            self.shape = (self.shape[0] - a.shape[0],) + self.shape[1:]
        else:
            assert a.shape == self.shape[1:]
            self.shape = (self.shape[0] - 1,) + self.shape[1:]
        assert self.shape[0] >= 0
            
        a.tofile(self.fd)

    def sync(self):
        """Write data dictionary.

        Write bool, int, float, complex and str data, shapes and
        dtypes for ndarrays and class names for other objects."""

        assert self.shape[0] == 0
        i = self.fd.tell()
        s = encode(self.data).encode()
        writeint(self.fd, len(s))
        self.fd.write(s)
        
        n = len(self.offsets)
        if self.nitems >= n:
            offsets = np.zeros(n * N1, np.int64)
            offsets[:n] = self.offsets
            self.itemoffsets = align(self.fd)
            offsets.tofile(self.fd)
            writeint(self.fd, self.itemoffsets, 24)
            self.offsets = offsets
            
        self.offsets[self.nitems] = i
        writeint(self.fd, i, self.itemoffsets + self.nitems * 8)
        self.nitems += 1
        writeint(self.fd, self.nitems, 16)
        self.fd.flush()
        self.fd.seek(0, 2)  # end of file
        self.data = {}
        
    def write(self, **kwargs):
        """Write data.

        Use::

            writer.write(n=7, s='abc', a=np.zeros(3), density=density).
        """
        
        for name, value in kwargs.items():
            if isinstance(value, (bool, int, float, complex,
                                  dict, list, tuple, str)):
                self.data[name] = value
            elif isinstance(value, np.ndarray):
                self.add_array(name, value.shape, value.dtype,
                               delayed_read=False)
                self.fill(value)
            else:
                self.data[name] = {'_type':
                                   value.__module__ + '.' +
                                   value.__class__.__name__}
                writer = Writer(self.fd, data=self.data[name])
                value.write(writer)
        
    def close(self):
        self.sync()
        self.fd.close()
        
        
def read_header(fd):
    fd.seek(0)
    assert fd.read(8) == b'IOASE...'
    version, nitems, itemoffsets = np.fromfile(fd, np.int64, 3)
    fd.seek(itemoffsets)
    offsets = np.fromfile(fd, np.int64, nitems)
    return version, nitems, itemoffsets, offsets

    
class Reader:
    def __init__(self, fd, item=0, data=None):
        """Create hierarchy of readers.

        Store data as attributes for easy access and to allow
        tab-completion."""
        
        assert np.little_endian

        if isinstance(fd, str):
            fd = open(fd, 'rb')
        
        self.fd = fd
        
        if data is None:
            self.version, self.nitems, self.itemoffsets, self.offsets = \
                read_header(fd)
            data = self._read_data(item)

        for name, value in data.items():
            if isinstance(value, dict) and '_type' in value:
                if value['_type'] == 'numpy.ndarray':
                    read_now = '_delayed' not in value
                    value = NDArrayReader(fd,
                                          value['shape'],
                                          np.dtype(value['dtype']),
                                          value['offset'])
                    if read_now:
                        value = value.read()
                else:
                    value = Reader(self.fd, data=value)
        
            data[name] = value
            
        self.data = data
        
    def __dir__(self):
        return self.data.keys()
    
    def __getattr__(self, attr):
        value = self.data[attr]
        if isinstance(value, NDArrayReader):
            return value.read()
        return value
        
    def proxy(self, name):
        value = self.data[name]
        assert isinstance(value, NDArrayReader)
        return value

    def __len__(self):
        return self.nitems
        
    def _read_data(self, item):
        self.fd.seek(self.offsets[item])
        size = np.fromfile(self.fd, np.int64, 1)[0]
        data = decode(self.fd.read(size).decode())
        return data
        
    def __getitem__(self, i):
        data = self._read_data(i)
        return Reader(self.fd, data=data)
        

read = Reader
write = Writer


class NDArrayReader:
    def __init__(self, fd, shape, dtype, offset):
        self.fd = fd
        self.shape = tuple(shape)
        self.dtype = dtype
        self.offset = offset
        
        self.ndim = len(self.shape)
        self.itemsize = dtype.itemsize
        self.size = np.prod(self.shape)
        self.nbytes = self.size * self.itemsize
        
    def __len__(self):
        return self.shape[0]
        
    def read(self):
        return self[:]
        
    def __getitem__(self, i):
        if isinstance(i, int):
            return self[i:i + 1][0]
        start, stop, step = i.indices(len(self))
        offset = self.offset + start * self.nbytes // len(self)
        self.fd.seek(offset)
        count = (stop - start) * self.size // len(self)
        a = np.fromfile(self.fd, self.dtype, count)
        a.shape = (-1,) + self.shape[1:]
        if step != 1:
            return a[::step].copy()
        return a

        
if __name__ == '__main__':
    class A:
        def write(self, writer):
            writer.write(x=np.ones((2, 3)))
    
        @staticmethod
        def read(reader):
            a = A()
            a.x = reader.x
            return a
        
    w = Writer('a.ioase')
    w.write(a=A(), y=9)
    w.write(s='abc')
    w.sync()
    w.write(s='abc2')
    w.sync()
    w.write(s='abc3', z=np.ones(7, int))
    w.close()
    print(w.data)
    
    r = Reader('a.ioase')
    print(r.y, r.s)
    print(A.read(r.a).x)
    print(r.a.x)
    print(r[1].s)
    print(r[2].s)
    print(r[2].z)
    
    w = Writer('a.ioase', 'a')
    print(w.nitems, w.offsets)
    w.write(d={'h': [1, 'asdf']})
    w.add_array('psi', (4, 3))
    w.fill(np.ones((1, 3)))
    w.fill(np.ones((1, 3)) * 2)
    w.fill(np.ones((2, 3)) * 3)
    w.close()
    print(Reader('a.ioase', 3).d)
    print(Reader('a.ioase')[2].z)
    print(Reader('a.ioase', 3).proxy('psi')[0:3])
