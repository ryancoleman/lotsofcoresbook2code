""" Test the integrate method of grid descriptor
    With this input, integrate uses dgemm for evaluation"""

from gpaw.grid_descriptor import GridDescriptor
from gpaw.test import equal
from gpaw.utilities.blas import gemm
from gpaw.mic.micblas import gemm as mic_gemm
from time import time
import pymic
import numpy as np
import sys

device = pymic.devices[0]
stream = device.get_default_stream()

gpts = int(sys.argv[1])
nbands = int(sys.argv[2])
repeats = 10

gd = GridDescriptor([gpts, gpts, gpts])

a = gd.empty(nbands)
b = gd.empty(nbands)

np.random.seed(10)
a[:] = np.random.random(a.shape)
c = np.random.random((nbands, nbands))

# warm-up
for i in range(3):
    gemm(1.0, a, c, 0.0, b)
# equal(np.sum(c), 3600.89536641, 1e-6)
t0 = time()
for i in range(repeats):
    gemm(1.0, a, c, 0.0, b)
t1 = time()
print "Check", np.sum(b), "Time", (t1 - t0) / repeats

a_mic = gd.empty(nbands, usemic=True)
b_mic = gd.empty(nbands, usemic=True)
np.random.seed(10)
a_mic.array[:] = np.random.random(a_mic.shape)
c_mic = device.associate(c)
a_mic.update_device()

# warm-up
for i in range(3):
    mic_gemm(1.0, a_mic, c_mic, 0.0, b_mic)
t0 = time()
# equal(np.sum(c), 3600.89536641, 1e-6)
for i in range(repeats):
    mic_gemm(1.0, a_mic, c_mic, 0.0, b_mic)
t1 = time()
b_mic.update_host()
print "Check", np.sum(b_mic.array), "Time", (t1 - t0) / repeats
