""" Test the integrate method of grid descriptor
    With this input, integrate uses dgemm for evaluation"""

from gpaw.grid_descriptor import GridDescriptor
from gpaw.test import equal
from time import time
import numpy as np
from gpaw.utilities.blas import gemm
from gpaw.mic.micblas import gemm as mic_gemm
#from gpaw.mic.micblas import r2k as mic_r2k
import pyMIC as mic
from time import time
import sys

device = mic.devices[0]

gpts = 48
nbands = 256
if len(sys.argv) > 1:
    nbands = int(sys.argv[1])

gd = GridDescriptor([gpts, gpts, gpts])

alpha = 1.0
beta = 0.0

repeats = 1

a = gd.empty(nbands)
# a = a.reshape((nbands, -1))
b = gd.empty(nbands)
# b = a.reshape((nbands, -1))
# c = np.zeros((a.shape[0], b.shape[0]))

np.random.seed(10)
a[:] = np.random.random(a.shape)
b[:] = np.random.random(b.shape)
c = np.zeros((len(a), len(b)))

t0 = time()
for r in range(repeats):
    gemm(alpha, a, b, beta, c, "c")
t1 = time()
print "Host time", t1 - t0

a_sum = np.sum(a)
b_sum = np.sum(b)
c_sum = np.sum(c)
print "Checks "
print "   sum(a)=" +str(a_sum)
print "   sum(b)=" +str(b_sum)
print "   sum(c)=" +str(c_sum)

a_mic = mic.offload_array(a.shape, dtype=float)
a_mic.fillfrom(a)
# a_mic.update_device()
# a_mic.array = a_mic.array.reshape(a_mic.shape[:-3] + (-1,))
a_mic.update_host()
a_mic_sum = np.sum(a_mic.array)
print "   sum(a_mic)=" + str(a_mic_sum)

b_mic = mic.offload_array(b.shape, dtype=float)
b_mic.fillfrom(b)
# b_mic.update_device()
# b_mic.array = b_mic.array.reshape(b_mic.shape[:-3] + (-1,))
b_mic.update_host()
b_mic_sum = np.sum(b_mic.array)
print "   sum(b_mic)=" + str(b_mic_sum)

# c_mic = offload_array(c.shape, dtype=float)
# c_mic.fill(0.0);
c_mic = device.associate(c)
c_mic.update_device()

t0 = time()
for r in range(repeats):
    mic_gemm(alpha, a_mic, b_mic, beta, c_mic, "c")
    #mic_r2k(alpha, a_mic, b_mic, beta, c_mic)
    c_mic.update_host()
    c[:] = c_mic.array[:]
t1 = time()
print "MIC time", t1 - t0
print "MIC checks"
c_mic_sum = np.sum(c_mic.array)
print "    sum(c_mic)=" + str(c_mic_sum)
print "    sum(c)=" + str(np.sum(c))
