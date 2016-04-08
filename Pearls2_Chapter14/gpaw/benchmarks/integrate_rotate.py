""" Test the integrate method of grid descriptor
    With this input, integrate uses dgemm for evaluation"""

from gpaw.grid_descriptor import GridDescriptor
from gpaw.test import equal
from time import time
import numpy as np
from gpaw.mpi import world, size, rank
from gpaw.utilities.blas import gemm
from gpaw.mic.micblas import gemm as mic_gemm
import pymic
import sys

device = pymic.devices[0]
stream = device.get_default_stream()

gpts = int(sys.argv[1])
nbands = int(sys.argv[2])
repeats = 10

gd = GridDescriptor([gpts, gpts, gpts], comm=world, parsize=size)

a = gd.empty(nbands)
b = gd.empty(nbands)

np.random.seed(10)
a[:] = np.random.random(a.shape)
b[:] = np.random.random(b.shape)
c = np.zeros((nbands, nbands))

# warm-up
for i in range(3):
    gd.integrate(a, b, hermitian=False, _transposed_result=c)
    gemm(1.0, a, c, 0.0, b)
# equal(np.sum(c), 3600.89536641, 1e-6)
t0 = time()
for i in range(repeats):
    gd.integrate(a, b, hermitian=False, _transposed_result=c)
    gemm(1.0, a, c, 0.0, b)
t1 = time()
if rank == 0:
    print "Check", np.sum(b), "Time", (t1 - t0) / repeats

a_mic = gd.empty(nbands, usemic=True)
b_mic = gd.empty(nbands, usemic=True)
c_mic = stream.bind(c)
np.random.seed(10)
a_mic.array[:] = np.random.random(a_mic.shape)
b_mic.array[:] = np.random.random(b_mic.shape)
# a_mic.update_device()
# b_mic.update_device()

# warm-up
for i in range(3):
    a_mic.update_device()
    b_mic.update_device()
    gd.integrate(a_mic, b_mic, hermitian=False, _transposed_result=c)
    c_mic.update_device()
    mic_gemm(1.0, a_mic, c_mic, 0.0, b_mic)
    b_mic.update_host()
t0 = time()
# equal(np.sum(c), 3600.89536641, 1e-6)
for i in range(repeats):
    a_mic.update_device()
    b_mic.update_device()
    gd.integrate(a_mic, b_mic, hermitian=False, _transposed_result=c)
    c_mic.update_device()
    mic_gemm(1.0, a_mic, c_mic, 0.0, b_mic)
    b_mic.update_host()
t1 = time()
if rank == 0:
    print "Check", np.sum(b_mic.array), "Time", (t1 - t0) / repeats
