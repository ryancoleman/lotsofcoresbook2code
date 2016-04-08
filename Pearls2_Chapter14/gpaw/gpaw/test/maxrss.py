from __future__ import print_function
import numpy as np
from gpaw.utilities.memory import maxrss, _VmB
from gpaw.test import equal

verbose = True

dtype = np.int32

def allocate_matrix(maxmem, msize, dtype):

    # square matrix in memory (in bytes): byte/bit * bit/int32 * size
    mmemory = (1./8 * 32 * msize**2)

    mno = int(max(1, maxmem * 2**20/mmemory)) # number of allocated matrices

    a = []
    for i in range(mno):
        a.append(np.ones((msize, msize), dtype=dtype))

    # must return and store the matrix to allocate memory
    return a, mmemory, mno

# max memory to be used for allocation matrices in MiB (2**20 B)
maxmems = [256] * 5 # exceed 1GiB to get resource.getrusage problem

msize = 256 # matrix dimensions

a = []

# initial memory
# Peak resident set size ("high water mark") (in bytes)
mem0 = maxrss()
mem = mem0

for n, maxmem in enumerate(maxmems):
    try:
        a1, mmemory, mno = allocate_matrix(maxmem, msize, dtype)
        # store the matrix
        a.append(a1)

        # memory used
        newmem = maxrss()
        memused = newmem - mem

        mset = 'Matrix set No ' + str(n) + ':'
        if verbose:
            print(mset, end=' ')
            print(str(mno) + ' matrices ' + str(dtype), end=' ')
            print('of size ' + str(msize) + '**2', end=' ')
            print('in memory of ' + str(int(mmemory*mno)) + ' bytes', end=' ')

            print('Total: ' + str(int(newmem - mem0)) + ' bytes')
 
        equal(float(memused)/(float(mmemory)*float(mno)), 1.0, 0.30, msg=mset)

    except MemoryError:
        if verbose:
            print(mset, end=' ')
            print('allocation failed') 
        pass
    mem = newmem

