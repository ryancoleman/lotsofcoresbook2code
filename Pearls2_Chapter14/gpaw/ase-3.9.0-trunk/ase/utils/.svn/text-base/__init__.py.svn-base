import os
import sys
import time
from math import sin, cos, radians, atan2, degrees

import numpy as np

from ase.data import chemical_symbols


class DevNull:
    def write(self, string):
        pass

    def flush(self):
        pass

    def seek(self, offset, whence=0):
        return 0

    def tell(self):
        return 0

    def close(self):
        pass
    
    def isatty(self):
        return False
        

devnull = DevNull()


def opencew(filename, world=None):
    """Create and open filename exclusively for writing.

    If master cpu gets exclusive write access to filename, a file
    descriptor is returned (a dummy file descriptor is returned on the
    slaves).  If the master cpu does not get write access, None is
    returned on all processors."""

    if world is None:
        from ase.parallel import world
        
    if world.rank == 0:
        try:
            fd = os.open(filename, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except OSError:
            ok = 0
        else:
            ok = 1
            fd = os.fdopen(fd, 'w')
    else:
        ok = 0
        fd = devnull

    # Syncronize:
    if world.sum(ok) == 0:
        return None
    else:
        return fd


class Lock:
    def __init__(self, name='lock', world=None):
        self.name = name
        
        if world is None:
            from ase.parallel import world
        self.world = world

    def acquire(self):
        while True:
            fd = opencew(self.name, self.world)
            if fd is not None:
                break
            time.sleep(1.0)
            
    def release(self):
        self.world.barrier()
        if self.world.rank == 0:
            os.remove(self.name)

    def __enter__(self):
        self.acquire()

    def __exit__(self, type, value, tb):
        self.release()


class OpenLock:
    def acquire(self):
        pass

    def release(self):
        pass

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        pass


def hill(numbers):
    """Convert list of atomic numbers to a chemical formula as a string.
    
    Elements are alphabetically ordered with C and H first."""
    
    d = {}
    for Z in numbers:
        symb = chemical_symbols[Z]
        d[symb] = d.get(symb, 0) + 1
    result = [(s, d.pop(s)) for s in 'CH' if s in d]
    result += [(s, d[s]) for s in sorted(d)]
    return ''.join('{0}{1}'.format(symbol, n) if n > 1 else symbol
                   for symbol, n in result)


def prnt(*args, **kwargs):
    """Python 3 style print function."""
    fd = kwargs.pop('file', sys.stdout)
    fd.write(
        kwargs.pop('sep', ' ').join(str(arg) for arg in args) +
        kwargs.pop('end', '\n'))
    if kwargs.pop('flush', False):
        fd.flush()
    if kwargs:
        raise TypeError('%r is an invalid keyword argument for this function' %
                        kwargs.keys()[0])


def gcd(a, b):
    """Greatest common divisor of a and b."""
    while a != 0:
        a, b = b % a, a
    return b


def rotate(rotations, rotation=np.identity(3)):
    """Convert string of format '50x,-10y,120z' to a rotation matrix.

    Note that the order of rotation matters, i.e. '50x,40z' is different
    from '40z,50x'.
    """

    if rotations == '':
        return rotation.copy()

    for i, a in [('xyz'.index(s[-1]), radians(float(s[:-1])))
                 for s in rotations.split(',')]:
        s = sin(a)
        c = cos(a)
        if i == 0:
            rotation = np.dot(rotation, [(1, 0, 0),
                                         (0, c, s),
                                         (0, -s, c)])
        elif i == 1:
            rotation = np.dot(rotation, [(c, 0, -s),
                                         (0, 1, 0),
                                         (s, 0, c)])
        else:
            rotation = np.dot(rotation, [(c, s, 0),
                                         (-s, c, 0),
                                         (0, 0, 1)])
    return rotation


def givens(a, b):
    """Solve the equation system::

      [ c s]   [a]   [r]
      [    ] . [ ] = [ ]
      [-s c]   [b]   [0]
    """
    sgn = lambda x: cmp(x, 0)
    if b == 0:
        c = sgn(a)
        s = 0
        r = abs(a)
    elif abs(b) >= abs(a):
        cot = a / b
        u = sgn(b) * (1 + cot**2)**0.5
        s = 1. / u
        c = s * cot
        r = b * u
    else:
        tan = b / a
        u = sgn(a) * (1 + tan**2)**0.5
        c = 1. / u
        s = c * tan
        r = a * u
    return c, s, r


def irotate(rotation, initial=np.identity(3)):
    """Determine x, y, z rotation angles from rotation matrix."""
    a = np.dot(initial, rotation)
    cx, sx, rx = givens(a[2, 2], a[1, 2])
    cy, sy, ry = givens(rx, a[0, 2])
    cz, sz, rz = givens(cx * a[1, 1] - sx * a[2, 1],
                        cy * a[0, 1] - sy * (sx * a[1, 1] + cx * a[2, 1]))
    x = degrees(atan2(sx, cx))
    y = degrees(atan2(-sy, cy))
    z = degrees(atan2(sz, cz))
    return x, y, z


def hsv2rgb(h, s, v):
    """http://en.wikipedia.org/wiki/HSL_and_HSV

    h (hue) in [0, 360[
    s (saturation) in [0, 1]
    v (value) in [0, 1]

    return rgb in range [0, 1]
    """
    if v == 0:
        return 0, 0, 0
    if s == 0:
        return v, v, v

    i, f = divmod(h / 60., 1)
    p = v * (1 - s)
    q = v * (1 - s * f)
    t = v * (1 - s * (1 - f))

    if i == 0:
        return v, t, p
    elif i == 1:
        return q, v, p
    elif i == 2:
        return p, v, t
    elif i == 3:
        return p, q, v
    elif i == 4:
        return t, p, v
    elif i == 5:
        return v, p, q
    else:
        raise RuntimeError('h must be in [0, 360]')


def hsv(array, s=.9, v=.9):
    array = (array + array.min()) * 359. / (array.max() - array.min())
    result = np.empty((len(array.flat), 3))
    for rgb, h in zip(result, array.flat):
        rgb[:] = hsv2rgb(h, s, v)
    return np.reshape(result, array.shape + (3,))

## This code does the same, but requires pylab
## def cmap(array, name='hsv'):
##     import pylab
##     a = (array + array.min()) / array.ptp()
##     rgba = getattr(pylab.cm, name)(a)
##     return rgba[:-1] # return rgb only (not alpha)

ON_POSIX = 'posix' in sys.builtin_module_names

try:
    from subprocess import Popen
except ImportError:
    from os import popen3
else:
    def popen3(cmd):
        from subprocess import PIPE
        p = Popen(cmd, shell=True, close_fds=ON_POSIX,
                  stdin=PIPE, stdout=PIPE, stderr=PIPE)
        return p.stdin, p.stdout, p.stderr
