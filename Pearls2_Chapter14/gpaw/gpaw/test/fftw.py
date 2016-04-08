from __future__ import print_function
import time

import numpy as np

import gpaw.fftw as fftw


def test(Plan, flags, input, output, sign):
    t0 = time.time()
    plan = Plan(input, output, sign, flags)
    t1 = time.time()
    t = 0.0
    for i in range(100):
        input[:] = 1.3
        t2 = time.time()
        plan.execute()
        t3 = time.time()
        t += t3 - t2
    return t1 - t0, t / 100


if __name__ == '__main__':
    a1 = fftw.empty((32, 28, 128), complex)
    a2 = fftw.empty((32, 28, 128), complex)
    b = fftw.empty((32, 28, 65), complex)
    c1 = b.view(dtype=float)[:, :, :128]
    c2 = fftw.empty((32, 28, 64), complex).view(dtype=float)
    for input, output, sign in [
        (a1, a1, -1),
        (a1, a2, -1),
        (b, c1, 1),
        (b, c2, 1),
        (c1, b, -1),
        (c2, b, -1)]:
        for Plan, flags in [(fftw.NumpyFFTPlan, 117),
                            (fftw.FFTWPlan, fftw.ESTIMATE),
                            (fftw.FFTWPlan, fftw.MEASURE),
                            (fftw.FFTWPlan, fftw.PATIENT),
                            (fftw.FFTWPlan, fftw.EXHAUSTIVE)]:
            tplan, tfft = test(Plan, flags, input, output, sign)
            print(('%-12s %3d %10.6f %10.6f' %
                  (Plan.__name__, flags, tplan, tfft)))
