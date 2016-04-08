# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""Utility functions and classes."""

import os
import re
import sys
from operator import mul
from math import sqrt, exp

import numpy as np
from numpy import linalg

import _gpaw
from gpaw import debug

from ase.utils import devnull

elementwise_multiply_add = _gpaw.elementwise_multiply_add
utilities_vdot = _gpaw.utilities_vdot
utilities_vdot_self = _gpaw.utilities_vdot_self


erf = np.vectorize(_gpaw.erf, (float,), 'Error function')
# XXX should we unify double and complex erf ???
cerf = np.vectorize(_gpaw.cerf, (complex,), 'Complex error function')

# Factorials:
_fact = [1, 1, 2, 6, 24, 120, 720, 5040, 40320,
         362880, 3628800, 39916800, 479001600]

def fact(x):
    if x in xrange(len(_fact)):
        return _fact[x]
    return reduce(mul, xrange(2, x+1), 1)

def ffact(a, b):
    """b!/a! where 0 <= a <= b"""
    assert a in xrange(b+1)
    return reduce(mul, xrange(a+1, b+1), 1)

# Code will crash for setups without any projectors.  Setups that have
# no projectors therefore receive a dummy projector as a hacky
# workaround.  The projector is assigned a certain, small size.  If
# the grid is so coarse that no point falls within the projector's range,
# there'll also be an error.  So this limits allowed grid spacings.
min_locfun_radius = 0.85 # Bohr
smallest_safe_grid_spacing = 2 * min_locfun_radius / np.sqrt(3) # ~0.52 Ang

def h2gpts(h, cell_cv, idiv=4):
    """Convert grid spacing to number of grid points divisible by idiv.

    Note that units of h and cell_cv must match!
    
    h: float
        Desired grid spacing in.
    cell_cv: 3x3 ndarray
        Unit cell.
    """

    L_c = (np.linalg.inv(cell_cv)**2).sum(0)**-0.5
    return np.maximum(idiv, (L_c / h / idiv + 0.5).astype(int) * idiv)


def gcd(a, b):
    """Return greatest common divisor of a and b, using the
    euclidian algorithm.
    """
    while b != 0:
        a, b = b, a % b
    return a


def is_contiguous(array, dtype=None):
    """Check for contiguity and type."""
    if dtype is None:
        return array.flags.c_contiguous
    else:
        return array.flags.c_contiguous and array.dtype == dtype


# Radial-grid Hartree solver:
#
#                       l
#             __  __   r
#     1      \   4||    <   * ^    ^
#   ------ =  )  ---- ---- Y (r)Y (r'),
#    _ _     /__ 2l+1  l+1  lm   lm
#   |r-r'|    lm      r
#                      >
# where
#
#   r = min(r, r')
#    <
#
# and
#
#   r = max(r, r')
#    >
#
def hartree(l, nrdr, r, vr):
    """Calculates radial Coulomb integral.

    The following integral is calculated::

                                   ^
                          n (r')Y (r')
              ^    / _     l     lm
      v (r)Y (r) = |dr' --------------,
       l    lm     /        _   _
                           |r - r'|

    where input and output arrays `nrdr` and `vr`::

              dr
      n (r) r --  and  v (r) r.
       l      dg        l
    """
    assert is_contiguous(nrdr, float)
    assert is_contiguous(r, float)
    assert is_contiguous(vr, float)
    assert nrdr.shape == vr.shape and len(vr.shape) == 1
    assert len(r.shape) == 1
    assert len(r) >= len(vr)
    return _gpaw.hartree(l, nrdr, r, vr)

    
def packed_index(i1, i2, ni):
    """Return a packed index"""
    if i1 > i2:
        return (i2 * (2 * ni - 1 - i2) // 2) + i1
    else:
        return (i1 * (2 * ni - 1 - i1) // 2) + i2


def unpacked_indices(p, ni):
    """Return unpacked indices corresponding to upper triangle"""
    assert 0 <= p < ni * (ni + 1) // 2
    i1 = int(ni + .5 - sqrt((ni - .5)**2 - 2 * (p - ni)))
    return i1, p - i1 * (2 * ni - 1 - i1) // 2


packing_conventions = """\n
In the code, the convention is that density matrices are constructed using
pack / unpack2, and anything that should be multiplied onto such, e.g.
corrections to the Hamiltonian, are constructed according to pack2 / unpack.
"""


def unpack(M):
    """Unpack 1D array to 2D, assuming a packing as in ``pack2``."""
    assert is_contiguous(M)
    n = int(sqrt(0.25 + 2.0 * len(M)))
    M2 = np.zeros((n, n), M.dtype.char)
    if M.dtype.char == complex:
        _gpaw.unpack_complex(M, M2)
    else:
        _gpaw.unpack(M, M2)
    return M2

    
def unpack2(M):
    """Unpack 1D array to 2D, assuming a packing as in ``pack``."""
    M2 = unpack(M)
    M2 *= 0.5 # divide all by 2
    M2.flat[0::len(M2) + 1] *= 2 # rescale diagonal to original size
    return M2

    
def pack(A):
    """Pack a 2D array to 1D, adding offdiagonal terms.
    
    The matrix::
    
           / a00 a01 a02 \ 
       A = | a10 a11 a12 |
           \ a20 a21 a22 /
                
    is transformed to the vector::
    
      (a00, a01 + a10, a02 + a20, a11, a12 + a21, a22)
    """
    assert A.ndim == 2
    assert A.shape[0] == A.shape[1]
    assert A.dtype in [float, complex]
    return _gpaw.pack(A)


def pack2(M2, tolerance=1e-10):
    """Pack a 2D array to 1D, averaging offdiagonal terms."""
    n = len(M2)
    M = np.zeros(n * (n + 1) // 2, M2.dtype.char)
    p = 0
    for r in range(n):
        M[p] = M2[r, r]
        p += 1
        for c in range(r + 1, n):
            M[p] = (M2[r, c] + np.conjugate(M2[c, r])) / 2. # note / 2.
            error = abs(M2[r, c] - np.conjugate(M2[c, r]))
            assert error < tolerance, 'Pack not symmetric by %s' % error + ' %'
            p += 1
    assert p == len(M)
    return M


for method in (unpack, unpack2, pack, pack2):
    method.__doc__ += packing_conventions


def element_from_packed(M, i, j):
    """Return a specific element from a packed array (by ``pack``)."""
    n = int(sqrt(2 * len(M) + .25)) 
    assert i < n and j < n
    p = packed_index(i, j, n)
    if i == j:
        return M[p]
    elif i > j:
        return .5 * M[p]
    else:
        return .5 * np.conjugate(M[p])
    

def logfile(name, rank=0):
    """Create file object from name.

    Use None for /dev/null and '-' for sys.stdout.  Ranks > 0 will
    get /dev/null."""

    if rank == 0:
        if name is None:
            fd = devnull
        elif name == '-':
            fd = sys.stdout
        elif isinstance(name, str):
            fd = open(name, 'w')
        else:
            fd = name
    else:
        fd = devnull
    return fd


def uncamelcase(name):
    """Convert a CamelCase name to a string of space-seperated words."""
    words = re.split('([A-Z]{1}[a-z]+)', name)
    return ' '.join([word for word in words if word != ''])


def divrl(a_g, l, r_g):
    """Return array divided by r to the l'th power."""
    b_g = a_g.copy()
    if l > 0:
        b_g[1:] /= r_g[1:]**l
        b1, b2 = b_g[1:3]
        r0, r1, r2 = r_g[0:3]
        b_g[0] = b2 + (b1 - b2) * (r0 - r2) / (r1 - r2)
    return b_g


def compiled_with_sl():
    return hasattr(_gpaw, 'new_blacs_context')


def load_balance(paw, atoms):
    try:
        paw.initialize(atoms)
    except SystemExit:
        pass
    spos_ac = paw.atoms.get_scaled_positions() % 1.0
    atoms_r = np.zeros(paw.wfs.world.size)
    rnk_a = paw.wfs.gd.get_ranks_from_positions(spos_ac)
    for rnk in rnk_a:
        atoms_r[rnk] += 1
    max_atoms = max(atoms_r)
    min_atoms = min(atoms_r)
    ave_atoms = atoms_r.sum()/paw.wfs.world.size
    stddev_atoms = sqrt((atoms_r**2).sum()/paw.wfs.world.size - ave_atoms**2)
    print("Information about load balancing")
    print("--------------------------------")
    print("Number of atoms:", len(spos_ac))
    print("Number of CPUs:", paw.wfs.world.size)
    print("Max. number of atoms/CPU:   ", max_atoms)
    print("Min. number of atoms/CPU:   ", min_atoms)
    print("Average number of atoms/CPU:", ave_atoms)
    print("    standard deviation:     %5.1f" % stddev_atoms)

if not debug:
    hartree = _gpaw.hartree
    pack = _gpaw.pack


def mlsqr(order, cutoff, coords_nc, N_c, beg_c, data_g, target_n):
    """Interpolate a point using moving least squares algorithm.

    Python wrapper for a c-function. See c/mlsqr.c.

    order:       Polynomial order (1 or 2)
    coords_nc:   List of scaled coordinates
    N_c:         Total number of grid points
    beg_c:       The start of grid points
    data_g:      3D-data to be interpolated
    target_n:    Output array
    """

    assert is_contiguous(coords_nc, float)
    assert is_contiguous(data_g, float)
    N_c = np.ascontiguousarray(N_c, float)
    beg_c = np.ascontiguousarray(beg_c, float)
    assert is_contiguous(target_n, float)

    return _gpaw.mlsqr(order, cutoff, coords_nc, N_c, beg_c, data_g, target_n)
    

def interpolate_mlsqr(dg_c, vt_g, order):
    """Interpolate a point using moving least squares algorithm.

    dg_c:    The grid point index (from (0,0,0) to (Bg_g - bg_g)
    vt_g:    The array to be interpolated
    order:   Polynomial order
    """

    # Define the weight function
    lsqr_weight = lambda r2: exp(-r2)

    # Define the polynomial basis
    if order == 1:
        b = lambda x: np.array([1, x[0], x[1], x[2]])
    elif order == 2:
        b = lambda x:  np.array([1, x[0], x[1], x[2],
                                 x[0] * x[1], x[1] * x[2], x[2] * x[0],
                                 x[0]**2, x[1]**2, x[2]**2])
    else:
        raise NotImplementedError

    def fill_X(x, y, z):
        result = None
        for i, j, k in zip(x.ravel(), y.ravel(), z.ravel()):
            r = b(np.array([i, j, k])) * lsqr_weight(
                np.sum((dg_c - np.array([i, j, k]))**2))
            if result is None:
                result = r
            else:
                result = np.vstack((result, r))
        return result

    def fill_w(x, y, z):
        result = []
        for i, j, k in zip(x.ravel(), y.ravel(), z.ravel()):
            weight = lsqr_weight(np.sum((dg_c - np.array([i, j, k]))**2))
            result.append(weight * vt_g[i][j][k])
        return np.array(result)
    
    X = np.fromfunction(fill_X, vt_g.shape)
    y = np.fromfunction(fill_w, vt_g.shape)

    X2 = np.dot(X.transpose(), X)
    y2 = np.dot(X.transpose(), y)
    c = linalg.solve(X2, y2)
    a = np.dot(b(dg_c), c)
    return a
