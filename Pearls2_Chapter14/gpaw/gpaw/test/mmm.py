"""Test BLAS matrix-matrix-multiplication interface."""
import numpy as np
from gpaw.utilities.blas import mmm


def op(o, m):
    if o == 'n':
        return m
    if o == 't':
        return m.T
    return m.T.conj()
    
    
def matrix(shape, dtype):
    if dtype == float:
        return np.random.random(shape)
    return np.random.random(shape) + 1j * np.random.random(shape)
    
        
for dtype in [float, complex]:
    a = matrix((2, 3), dtype)
    for opa in 'ntc':
        A = op(opa, a)
        B = matrix((A.shape[1], 4), dtype)
        for opb in 'ntc':
            b = op(opb, B).copy()
            C = np.dot(A, B)
            mmm(1, a, opa, b, opb, -1, C)
            assert abs(C).max() < 1e-14
