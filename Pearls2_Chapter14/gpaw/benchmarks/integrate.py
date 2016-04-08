""" Test the integrate method of grid descriptor
    With this input, integrate uses dgemm for evaluation"""

from gpaw.grid_descriptor import GridDescriptor
from gpaw.test import equal
from time import time
import numpy as np
import sys

gpts = int(sys.argv[1])
nbands = int(sys.argv[2])
repeats = 10

gd = GridDescriptor([gpts, gpts, gpts])

a = gd.empty(nbands)
b = gd.empty(nbands)

np.random.seed(10)
a[:] = np.random.random(a.shape)
b[:] = np.random.random(b.shape)

# warm-up
for i in range(3):
    c = gd.integrate(a, b)
# equal(np.sum(c), 3600.89536641, 1e-6)
t0 = time()
for i in range(repeats):
    c = gd.integrate(a, b)
t1 = time()
print "Check", np.sum(c), "Time", (t1 - t0) / repeats

a_mic = gd.empty(nbands, usemic=True)
b_mic = gd.empty(nbands, usemic=True)
np.random.seed(10)
a_mic.array[:] = np.random.random(a_mic.shape)
b_mic.array[:] = np.random.random(b_mic.shape)
a_mic.update_device()
b_mic.update_device()

# warm-up
for i in range(3):
    c = gd.integrate(a_mic, b_mic)
t0 = time()
# equal(np.sum(c), 3600.89536641, 1e-6)
for i in range(repeats):
    c = gd.integrate(a_mic, b_mic)
t1 = time()
print "Check", np.sum(c), "Time", (t1 - t0) / repeats
