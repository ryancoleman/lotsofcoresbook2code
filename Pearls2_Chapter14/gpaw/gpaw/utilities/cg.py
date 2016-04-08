import numpy as np

def CG(A, X, B, maxiter=20, tolerance=1.0e-10, verbose=False):
    """Solve X*A=B using conjugate gradient method.

    ``X`` and ``B`` are ``ndarrays```of shape ``(m, nx, ny, nz)``
    coresponding to matrices of size ``m*n`` (``n=nx*ny*nz``) and
    ``A`` is a callable representing an ``n*n`` matrix::

      A(X, Y)

    will store ``X*A`` in the output array ``Y``.
    
    On return ``X`` will be the solution to ``X*A=B`` within
    ``tolerance``."""

    m = len(X)
    shape = (m, 1, 1, 1)
    R = np.empty(X.shape, X.dtype.char)
    Q = np.empty(X.shape, X.dtype.char)
    A(X, R)
    R -= B
    P = R.copy()
    c1 = A.sum(np.reshape([abs(np.vdot(r, r)) for r in R], shape))
    for i in range(maxiter):
        error = sum(c1.ravel())
        if verbose:
            print('CG-%d: %e' % (i, error))
        if error < tolerance:
            return i, error
        A(P, Q)
        #alpha = c1 / reshape([vdot(p, q) for p, q in zip(P, Q)], shape)
        alpha = c1 / A.sum(np.reshape([np.vdot(q,p)
                                       for p, q in zip(P, Q)], shape))
        X -= alpha * P
        R -= alpha * Q
        c0 = c1
        c1 = A.sum(np.reshape([abs(np.vdot(r, r)) for r in R], shape))
        beta = c1 / c0
        P *= beta
        P += R
        
    raise ArithmeticError('Did not converge!')
