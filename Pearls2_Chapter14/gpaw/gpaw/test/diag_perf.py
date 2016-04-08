from __future__ import print_function
import numpy as np
from time import time

from gpaw.utilities.lapack import diagonalize, diagonalize_mr3

seed = 43
gen = np.random.RandomState(seed)

def main(i,seed=42,dtype=float):
    if (dtype==complex):
        epsilon = 0.1j
    else:
        epsilon = 0.1
    x = i + 1
    N = x*100
    print("N =",N)
    H0 = np.zeros((N,N),dtype=dtype) + gen.rand(*(N,N))
    H1 = H0 + epsilon*np.tri(N,N, k=-1)
    W0 = np.zeros((N))
    Z0 = np.zeros_like(H0)
    t0 = time()
    diagonalize(H1, W0)
    t1 = time() - t0
    print("diagonalize", t1)
    t2 = time() 
    diagonalize_mr3(H1, W0, Z0)
    t3 = time() - t2
    print("diagonalize_mr3",t3)
    # diagonalize_mr3 must be faster than diagonalize 
    assert(t3 < t1)

if __name__ in ['__main__', '__builtin__']:
    for i in range(8): # Check matrix sizes only up to 800
        main(i,dtype=float)
        main(i,dtype=complex)
