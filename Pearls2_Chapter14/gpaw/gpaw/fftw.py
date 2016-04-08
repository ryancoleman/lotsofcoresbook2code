"""Python wrapper for FFTW3 library."""

import os

import numpy as np

if 'GPAW_FFTWSO' in os.environ:
    fftwlibnames = [os.environ['GPAW_FFTWSO']]
    if '' in fftwlibnames:
        fftwlibnames.remove('')  # GPAW_FFTWSO='' for numpy fallback!
else:
    fftwlibnames = ['libmkl_rt.so', 'libmkl_intel_lp64.so', 'libfftw3.so']

lib = None
for libname in fftwlibnames:
    try:
        import ctypes
        lib = ctypes.CDLL(libname)
        break
    except (ImportError, OSError):
        pass

ESTIMATE = 64
MEASURE = 0
PATIENT = 32
EXHAUSTIVE = 8


def check_fft_size(n):
    """Check if n is an efficient fft size.

    Efficient means that n can be factored into small primes (2, 3, 5, 7)."""

    if n == 1:
        return True
    for x in [2, 3, 5, 7]:
        if n % x == 0:
            return check_fft_size(n // x)
    return False


def get_efficient_fft_size(N, n=1):
    """Return smallest efficient fft size.
    
    Must be greater than or equal to N and divisible by n.
    """
    N = -(-N // n) * n
    while not check_fft_size(N):
        N += n
    return N


class FFTWPlan:
    """FFTW3 3d transform."""
    def __init__(self, in_R, out_R, sign, flags=MEASURE):
        if in_R.dtype == float:
            assert sign == -1
            n0, n1, n2 = in_R.shape
            self.plan = lib.fftw_plan_dft_r2c_3d(n0, n1, n2,
                                                 in_R, out_R, flags)
        elif out_R.dtype == float:
            assert sign == 1
            n0, n1, n2 = out_R.shape
            self.plan = lib.fftw_plan_dft_c2r_3d(n0, n1, n2,
                                                 in_R, out_R, flags)
        else:
            n0, n1, n2 = in_R.shape
            self.plan = lib.fftw_plan_dft_3d(n0, n1, n2,
                                             in_R, out_R, sign, flags)

    def execute(self):
        lib.fftw_execute(self.plan)

    def __del__(self, lib=lib):
        lib.fftw_destroy_plan(self.plan)


class NumpyFFTPlan:
    """Numpy fallback."""
    def __init__(self, in_R, out_R, sign, flags=None):
        self.in_R = in_R
        self.out_R = out_R
        self.sign = sign

    def execute(self):
        if self.in_R.dtype == float:
            self.out_R[:] = np.fft.rfftn(self.in_R)
        elif self.out_R.dtype == float:
            self.out_R[:] = np.fft.irfftn(self.in_R, self.out_R.shape)
            self.out_R *= self.out_R.size
        elif self.sign == 1:
            self.out_R[:] = np.fft.ifftn(self.in_R, self.out_R.shape)
            self.out_R *= self.out_R.size
        else:
            self.out_R[:] = np.fft.fftn(self.in_R)


def empty(shape, dtype=float):
    """numpy.empty() equivalent with 16 byte allignment."""
    assert dtype == complex
    N = np.prod(shape)
    a = np.empty(2 * N + 1)
    offset = (a.ctypes.data % 16) // 8
    a = a[offset:2 * N + offset].view(complex)
    a.shape = shape
    return a


if lib is None:
    FFTPlan = NumpyFFTPlan
else:
    FFTPlan = FFTWPlan
    lib.fftw_plan_dft_3d.argtypes = [
        ctypes.c_int, ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags='C_CONTIGUOUS'),
        ctypes.c_int, ctypes.c_uint]
    lib.fftw_plan_dft_r2c_3d.argtypes = [
        ctypes.c_int, ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=float, ndim=3),
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags='C_CONTIGUOUS'),
        ctypes.c_uint]
    lib.fftw_plan_dft_c2r_3d.argtypes = [
        ctypes.c_int, ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=float, ndim=3),
        ctypes.c_uint]
    try:
        lib.fftw_plan_with_nthreads.argtypes = [ctypes.c_int]
        assert lib.fftw_init_threads()
    except AttributeError:
        pass
