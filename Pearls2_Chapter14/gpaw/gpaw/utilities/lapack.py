# Copyright (C) 2003  CAMP
# Copyright (C) 2010  Argonne National Laboratory
# Please see the accompanying LICENSE file for further information.

"""
Python wrapper functions for the ``C`` package:
Linear Algebra PACKage (LAPACK)
"""

import numpy as np

from gpaw import debug
import _gpaw
from gpaw.utilities.tools import tri2full
from gpaw.utilities.blas import gemm


def diagonalize(a, w):
    """Diagonalize a symmetric/hermitian matrix.

    Uses dsyevd/zheevd to diagonalize symmetric/hermitian matrix
    `a`. The eigenvectors are returned in the rows of `a`, and the
    eigenvalues in `w` in ascending order. Only the lower triangle of
    `a` is considered."""

    assert a.flags.contiguous
    assert w.flags.contiguous
    assert a.dtype in [float, complex]
    assert w.dtype == float
    n = len(a)
    assert a.shape == (n, n)
    assert w.shape == (n,)

    info = _gpaw.diagonalize(a, w)
    if info != 0:
        raise RuntimeError('diagonalize error: %d' % info)


def diagonalize_mr3(a, w, z):
    """Diagonalize a symmetric/hermitian matrix.

    Uses dsyevr/zheevr to diagonalize symmetric/hermitian matrix
    `a`. The eigenvectors are returned in the rows of `z`, and the
    eigenvalues in `w` in ascending order. Only the lower triangle of
    `a` is considered."""

    assert a.flags.contiguous
    assert w.flags.contiguous
    assert z.flags.contiguous
    assert a.dtype in [float, complex]
    assert w.dtype == float
    assert z.dtype == a.dtype
    n = len(a)
    assert a.shape == (n, n)
    assert w.shape == (n,)
    assert z.shape == (n, n)
    info = _gpaw.diagonalize_mr3(a, w, z)
    if info != 0:
        raise RuntimeError('diagonalize_mr3 error: %d' % info)


def general_diagonalize(a, w, b, iu=None):
    """Diagonalize a generalized symmetric/hermitian matrix
    
    A * x = (lambda) * B * x,

    where `lambda` is the eigenvalue and `A` and `B` are the 
    matrices corresponding to `a` and `b`, respectively.

    If `iu` is `None`:
    Uses dsygvd/zhegvd to diagonalize symmetric/hermitian matrix
    `a`. The eigenvectors are returned in the rows of `a`, and the
    eigenvalues in `w` in ascending order. Only the lower triangle of
    `a` is considered.

    If `iu` is not `None`:
    Uses dsygvx/zhegvx to find the eigenvalues of 1 through `iu`.
    Stores the eigenvectors in `z` and the eigenvalues in `w` in
    ascending order.
    """

    assert a.flags.contiguous
    assert w.flags.contiguous
    assert a.dtype in [float, complex]
    assert w.dtype == float
    n = len(a)
    assert a.shape == (n, n)
    assert w.shape == (n,)
    assert b.flags.contiguous
    assert b.dtype == a.dtype
    assert b.shape == a.shape

    if iu is not None:
        z = np.zeros((n, n), dtype=a.dtype)
        assert z.flags.contiguous

    w[:1] = 42

    if iu is None:
        info = _gpaw.general_diagonalize(a, w, b)
    else:
        info = _gpaw.general_diagonalize(a, w, b, z, iu)
        a[:] = z

    if info != 0 or n > 0 and w[0] == 42:
        raise RuntimeError('general_diagonalize error: %d' % info)


def inverse_cholesky(a):
    """Calculate the inverse of the Cholesky decomposition of
    a symmetric/hermitian positive definite matrix `a`.

    Uses dpotrf/zpotrf to calculate the decomposition and then
    dtrtri/ztrtri for the inversion"""

    assert a.flags.contiguous
    assert a.dtype in [float, complex]
    n = len(a)
    assert a.shape == (n, n)

    info = _gpaw.inverse_cholesky(a)
    if info != 0:
        raise RuntimeError('inverse_cholesky error: %d' % info)


def inverse_general(a):
    assert a.dtype in [float, complex]
    n = len(a)
    assert a.shape == (n, n)
    info = _gpaw.inverse_general(a)
    if info != 0:
        raise RuntimeError('inverse_general error: %d' % info)


def inverse_symmetric(a):
    assert a.dtype in [float, complex]
    n = len(a)
    assert a.shape == (n, n)
    info = _gpaw.inverse_symmetric(a)
    tri2full(a, 'L', 'symm')
    if info != 0:
        raise RuntimeError('inverse_symmetric: %d' % info)


def right_eigenvectors(a, w, v):
    """Get right eigenvectors and eigenvalues from a square matrix
    using LAPACK dgeev.

    The right eigenvector corresponding to eigenvalue w[i] is v[i]."""

    assert a.flags.contiguous
    assert w.flags.contiguous
    assert v.flags.contiguous
    assert a.dtype == float
    assert w.dtype == float
    assert v.dtype == float
    n = len(a)
    assert a.shape == (n, n)
    assert w.shape == (n,)
    assert w.shape == (n, n)
    return _gpaw.right_eigenvectors(a, w, v)


def pm(M):
    """print a matrix or a vector in mathematica style"""
    string = ''
    s = M.shape
    if len(s) > 1:
        n, m = s
        string += '{'
        for i in range(n):
            string += '{'
            for j in range(m):
                string += str(M[i, j])
                if j == m - 1:
                    string += '}'
                else:
                    string += ','
            if i == n - 1:
                string += '}'
            else:
                string += ','
    else:
        n = s[0]
        string += '{'
        for i in range(n):
            string += str(M[i])
            if i == n - 1:
                string += '}'
            else:
                string += ','
    return string


def sqrt_matrix(a, preserve=False):
    """Get the sqrt of a symmetric matrix a (diagonalize is used).
    The matrix is kept if preserve=True, a=sqrt(a) otherwise."""
    n = len(a)
    if debug:
        assert a.flags.contiguous
        assert a.dtype == float
        assert a.shape == (n, n)
    if preserve:
        b = a.copy()
    else:
        b = a

    # diagonalize to get the form b = Z * D * Z^T
    # where D is diagonal
    D = np.empty((n,))
    diagonalize(b, D)
    ZT = b.copy()
    Z = np.transpose(b)

    # c = Z * sqrt(D)
    c = Z * np.sqrt(D)

    # sqrt(b) = c * Z^T
    gemm(1., ZT, np.ascontiguousarray(c), 0., b)
    return b
