from __future__ import print_function

import time
import numpy as np
from gpaw.utilities.gauss import gaussian_wave

sigma = 2.4
C = 3
G = 50

t = time.time()
r_cG = np.random.normal(size=C*G**3).reshape((C,G,G,G))
r0_c = np.random.normal(size=C)
k_c = np.random.normal(size=C)
A = np.random.uniform()*np.exp(1j*np.random.uniform(0,2*np.pi))
print('Allocation: %8.5f s' % (time.time()-t))

# -------------------------------------------------------------------

# Test case for real-part of gamma-point wave with normalized amplitude
_gaussRGN = lambda r_cG, r0_c, sigma: 1/(sigma*np.pi**0.5)**1.5 \
    * np.exp(-np.sum((r_cG-r0_c[:,np.newaxis,np.newaxis,np.newaxis])**2, axis=0)/(2*sigma**2))

t = time.time()
gs0_G = _gaussRGN(r_cG, r0_c, sigma)
print('_gaussRGN: %8.5f s' % (time.time()-t))

t = time.time()
gs1_G = gaussian_wave(r_cG, r0_c, sigma)
print('+gaussRGN: %8.5f s' % (time.time()-t))

assert np.abs(gs0_G-gs1_G).max() < 1e-12, 'Max error %g' % np.abs(gs0_G-gs1_G).max()
del gs0_G, gs1_G

# Test case for real-part of gamma-point wave with complex amplitude
_gaussRGA = lambda r_cG, r0_c, sigma, A: np.real(A) \
    * np.exp(-np.sum((r_cG-r0_c[:,np.newaxis,np.newaxis,np.newaxis])**2, axis=0)/(2*sigma**2))

t = time.time()
gs0_G = _gaussRGA(r_cG, r0_c, sigma, A)
print('_gaussRGA: %8.5f s' % (time.time()-t))

t = time.time()
gs1_G = gaussian_wave(r_cG, r0_c, sigma, None, A)
print('+gaussRGA: %8.5f s' % (time.time()-t))

assert np.abs(gs0_G-gs1_G).max() < 1e-12, 'Max error %g' % np.abs(gs0_G-gs1_G).max()
del gs0_G, gs1_G

# Test case for real-part of kpoint-point wave with normalized amplitude
_gaussRKN = lambda r_cG, r0_c, sigma, k_c: 1/(sigma*np.pi**0.5)**1.5 \
    * np.exp(-np.sum((r_cG-r0_c[:,np.newaxis,np.newaxis,np.newaxis])**2, axis=0)/(2*sigma**2)) \
    * np.cos(np.sum(r_cG*k_c[:,np.newaxis,np.newaxis,np.newaxis], axis=0))

t = time.time()
gs0_G = _gaussRKN(r_cG, r0_c, sigma, k_c)
print('_gaussRKN: %8.5f s' % (time.time()-t))

t = time.time()
gs1_G = gaussian_wave(r_cG, r0_c, sigma, k_c)
print('+gaussRKN: %8.5f s' % (time.time()-t))

assert np.abs(gs0_G-gs1_G).max() < 1e-12, 'Max error %g' % np.abs(gs0_G-gs1_G).max()
del gs0_G, gs1_G

# Test case for real-part of kpoint-point wave with complex amplitude
_gaussRKA = lambda r_cG, r0_c, sigma, k_c, A: \
    np.exp(-np.sum((r_cG-r0_c[:,np.newaxis,np.newaxis,np.newaxis])**2, axis=0)/(2*sigma**2)) \
    * np.real(A*np.exp(1j*np.sum(r_cG*k_c[:,np.newaxis,np.newaxis,np.newaxis], axis=0)))

t = time.time()
gs0_G = _gaussRKA(r_cG, r0_c, sigma, k_c, A)
print('_gaussRKA: %8.5f s' % (time.time()-t))

t = time.time()
gs1_G = gaussian_wave(r_cG, r0_c, sigma, k_c, A)
print('+gaussRKA: %8.5f s' % (time.time()-t))

assert np.abs(gs0_G-gs1_G).max() < 1e-12, 'Max error %g' % np.abs(gs0_G-gs1_G).max()
del gs0_G, gs1_G

# -------------------------------------------------------------------

# Test case for complex case of gamma-point wave with normalized amplitude
_gaussCGN = lambda r_cG, r0_c, sigma: (1+0j)/(sigma*np.pi**0.5)**1.5 \
    * np.exp(-np.sum((r_cG-r0_c[:,np.newaxis,np.newaxis,np.newaxis])**2, axis=0)/(2*sigma**2))

t = time.time()
gs0_G = _gaussCGN(r_cG, r0_c, sigma)
print('_gaussCGN: %8.5f s' % (time.time()-t))

t = time.time()
gs1_G = gaussian_wave(r_cG, r0_c, sigma, dtype=complex)
print('+gaussCGN: %8.5f s' % (time.time()-t))

assert np.abs(gs0_G-gs1_G).max() < 1e-12, 'Max error %g' % np.abs(gs0_G-gs1_G).max()
del gs0_G, gs1_G

# Test case for complex case of gamma-point wave with complex amplitude
_gaussCGA = lambda r_cG, r0_c, sigma, A: A \
    * np.exp(-np.sum((r_cG-r0_c[:,np.newaxis,np.newaxis,np.newaxis])**2, axis=0)/(2*sigma**2))

t = time.time()
gs0_G = _gaussCGA(r_cG, r0_c, sigma, A)
print('_gaussCGA: %8.5f s' % (time.time()-t))

t = time.time()
gs1_G = gaussian_wave(r_cG, r0_c, sigma, None, A, dtype=complex)
print('+gaussCGA: %8.5f s' % (time.time()-t))

assert np.abs(gs0_G-gs1_G).max() < 1e-12, 'Max error %g' % np.abs(gs0_G-gs1_G).max()
del gs0_G, gs1_G

# Test case for complex case of kpoint-point wave with normalized amplitude
_gaussCKN = lambda r_cG, r0_c, sigma, k_c: (1+0j)/(sigma*np.pi**0.5)**1.5 \
    * np.exp(-np.sum((r_cG-r0_c[:,np.newaxis,np.newaxis,np.newaxis])**2, axis=0)/(2*sigma**2)) \
    * np.exp(1j*np.sum(r_cG*k_c[:,np.newaxis,np.newaxis,np.newaxis], axis=0))

t = time.time()
gs0_G = _gaussCKN(r_cG, r0_c, sigma, k_c)
print('_gaussCKN: %8.5f s' % (time.time()-t))

t = time.time()
gs1_G = gaussian_wave(r_cG, r0_c, sigma, k_c, dtype=complex)
print('+gaussCKN: %8.5f s' % (time.time()-t))

assert np.abs(gs0_G-gs1_G).max() < 1e-12, 'Max error %g' % np.abs(gs0_G-gs1_G).max()
del gs0_G, gs1_G

# Test case for complex case of kpoint-point wave with complex amplitude
_gaussCKA = lambda r_cG, r0_c, sigma, k_c, A: A \
    * np.exp(-np.sum((r_cG-r0_c[:,np.newaxis,np.newaxis,np.newaxis])**2, axis=0)/(2*sigma**2)) \
    * np.exp(1j*np.sum(r_cG*k_c[:,np.newaxis,np.newaxis,np.newaxis], axis=0))

t = time.time()
gs0_G = _gaussCKA(r_cG, r0_c, sigma, k_c, A)
print('_gaussCKA: %8.5f s' % (time.time()-t))

t = time.time()
gs1_G = gaussian_wave(r_cG, r0_c, sigma, k_c, A, dtype=complex)
print('+gaussCKA: %8.5f s' % (time.time()-t))

assert np.abs(gs0_G-gs1_G).max() < 1e-12, 'Max error %g' % np.abs(gs0_G-gs1_G).max()
del gs0_G, gs1_G

