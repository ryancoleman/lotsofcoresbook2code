""" Test the Jacobi relaxation"""

from gpaw.grid_descriptor import GridDescriptor
from gpaw.fd_operators import Laplace
from gpaw.test import equal
import numpy as np
import time
import timeit
import os

import pyMIC as mic

offload_enabled = os.environ.get("GPAW_OFFLOAD");

device = mic.devices[0]

gpts = 64
gd = GridDescriptor([gpts, gpts, gpts])

a = gd.empty() # source
b = gd.empty() # result

np.random.seed(10)
a[:] = np.random.random(a.shape)
b[:] = np.zeros(b.shape)

op = Laplace(gd, 1.0, 3)

if offload_enabled:
    offl_a = device.bind(a)
    offl_b = device.bind(b)

print "--------------------------------------------------------------"
print "Starting relaxation step"
relax_start = time.time()
if offload_enabled:
    op.relax(2, offl_b.array, offl_a.array, 4, 2.0 / 3.0)
else:
    op.relax(2, b, a, 4, 2.0 / 3.0)
relax_end = time.time()
print "--------------------------------------------------------------"

print "time needed: " + str(relax_end-relax_start)

s = np.sum(b)
# print "Checksum: {0}".format(s)
if gpts == 64:
    equal(s, -10.4429888082, 1e-5)
if gpts == 128:
    equal(s, -20.9047491129, 1e-5)
    
#if __name__ == '__main__':
#    print timeit.repeat('op.relax(2, b, a, 4, 2.0 / 3.0)', setup='from  __main__ import op, a, b', number=100, repeat=5)

