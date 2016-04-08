import numpy as np
from numpy.fft import fft, ifft

"""
Methods for perfoming the Hilbert transform of a function::

                     +oo
            1       /     f(x) 
  H[f](y) = -- p.v. | dx -----
            pi      /    x - y  
                     -oo
"""

def hilbert_kernel_simple(n):
    """Construct Hilbert kernel with n grid points.
    
    This is just the discrete Fourier transform of 1 / x.
    """
    ker = np.zeros(n, dtype=complex)
    ker[1: n / 2] = 1.j
    ker[n / 2 + 1:] = -1.j
    return ker

def hilbert_kernel_interpolate(n):
    """Construct Hilbert kernel with n grid points.
    
    This is just the discrete Hilbert transform of the linear
    interpolation kernel `L(s) = (1 - |s|) Heaviside(1 - |s|)`.
    """
    # middle grid point
    mid = (n + 1) / 2

    # Make auxiliary array
    aux = np.arange(mid + 1, dtype=float)
    np.multiply(aux[1:], np.log(aux[1:]), aux[1:])

    # Make kernel
    ker = np.zeros(n, float)
    ker[1: mid] = aux[2:] - 2 * aux[1:-1] + aux[:-2]
    ker[-1: -mid: -1] = -ker[1: mid]
    
    return -fft(ker) / np.pi

def hilbert(f, ker=None, nfft=None, axis=0,
            kerneltype='interpolate', translate=0):
    """Perform Hilbert transform *f* along specified *axis*. The transform is
       made as a convolution of *f* with the Hilbert kernel.
       
       *ker* is the Hilbert kernel, which will be calculated if set to None.
       
       *nfft* is the number of grid points used in the fft. If *nfft* is larger
       than *f* along the transform axis, *f* will be zero-padded to make up
       the difference. If *nfft* is smaller, the first *nfft* elements of *f*
       will be used. *nfft* defaults to two times the the length of *f* along
       *axis*.
       
       *kerneltype* specifies the kerneltype, and can be one of 'interpolate'
       or 'simple'.
       
       *translate* is the number of grid points *f* should be shifted. This can
       be a non-integer amount. translate=10 means that f(x + 10) is
       transformed instead of f(x).
    """
    # Number of transform grid points
    n = f.shape[axis]

    # Number of grid points in fft
    if nfft is None: nfft = 2 * n

    # Generate new kernel if needed
    if ker is None: ker = eval('hilbert_kernel_' + kerneltype)(nfft)

    # Reshape kernel
    ker_shape = [1,] * len(f.shape)
    ker_shape[axis] = nfft
    ker.shape = tuple(ker_shape)

    # Construct translation operator
    if translate == 0:
        trans = 1
    else:
        trans = (np.arange(nfft) + (nfft - 1) / 2) % nfft - (nfft - 1) / 2
#          nfft = 8 ->       [ 0  1  2  3  4 -3 -2 -1]
#       trans = (np.arange(nfft) + nfft / 2) % nfft - nfft / 2
#          nfft = 8 ->       [ 0  1  2  3 -4 -3 -2 -1]
        trans = np.exp(1.j * translate * 2 * np.pi / nfft * trans)
        trans.shape = ker.shape

    # Make convolution of f and kernel
    hil = ifft(fft(f, n=nfft, axis=axis) * ker * trans, axis=axis)
    return hil[0:n]
    
def analytic_transforms(x):
    """Returns a list of functions and their Hilbert transforms"""
    func_l = [np.sin(x),
              np.cos(x),
              np.where(x == 0, 1, np.sin(x) / x),
              1 / (1 + x**2)]

    hfunc_l = [np.cos(x),
               -np.sin(x),
               np.where(x == 0, 0, (np.cos(x) - 1) / x),
               -x / (1 + x**2)]

    return func_l, hfunc_l
