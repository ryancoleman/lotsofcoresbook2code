# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

import numpy as np

from gpaw.utilities import erf


def I(R, a, b, alpha, beta):
    """Calculate integral and derivatives wrt. positions of Gaussian product.

    ::

                /                              2             2
               |  a0+b0  a1+b1  a2+b2  -alpha r   -beta (r-R)
       value = | x      y      z      e          e             dr,
               |
              /

    Returns the tuple (value, d[value]/dx, d[value]/dy, d[value]/dz).
    """
    result = np.zeros(4)
    R = np.array(R)
    result[0] = I1(R, a, b, alpha, beta)
    a = np.array(a)
    for i in range(3):
        a[i] += 1
        result[1 + i] = 2 * alpha * I1(R, tuple(a), b, alpha, beta)
        a[i] -= 2
        if a[i] >= 0:
            result[1 + i] -= (a[i] + 1) * I1(R, tuple(a), b, alpha, beta)
        a[i] += 1
    return result


def I1(R, ap1, b, alpha, beta, m=0):
    if ap1 == (0, 0, 0):
        if b != (0, 0, 0):
            return I1(-R, b, ap1, beta, alpha, m)
        else:
            f = 2 * np.sqrt(np.pi**5 / (alpha + beta)) / (alpha * beta)
            if np.sometrue(R):
                T = alpha * beta / (alpha + beta) * np.dot(R, R)
                f1 = f * erf(T**0.5) * (np.pi / T)**0.5
                if m == 0:
                    return 0.5 * f1
                f2 = f * np.exp(-T) / T**m
                if m == 1:
                    return 0.25 * f1 / T - 0.5 * f2
                if m == 2:
                    return 0.375 * f1 / T**2 - 0.5 * f2 * (1.5 + T)
                if m == 3:
                    return 0.9375 * f1 / T**3 - 0.25 * f2 * (7.5 +
                                                             T * (5 + 2 * T))
                if m == 4:
                    return 3.28125 * f1 / T**4 - 0.125 * f2 * \
                           (52.5 + T * (35 + 2 * T * (7 + 2 * T)))
                if m == 5:
                    return 14.7656 * f1 / T**5 - 0.03125 * f2 * \
                           (945 + T * (630 + T * (252 + T * (72 + T * 16))))
                if m == 6:
                    return 81.2109 * f1 / T**6 - 0.015625 * f2 * \
                           (10395 + T *
                            (6930 + T *
                             (2772 + T * (792 + T * (176 + T * 32)))))
                else:
                    raise NotImplementedError
                
            return f / (1 + 2 * m)
    for i in range(3):
        if ap1[i] > 0:
            break
    a = ap1[:i] + (ap1[i] - 1,) + ap1[i + 1:]
    result = beta / (alpha + beta) * R[i] * I1(R, a, b, alpha, beta, m + 1)
    if a[i] > 0:
        am1 = a[:i] + (a[i] - 1,) + a[i + 1:]
        result += a[i] / (2 * alpha) * (I1(R, am1, b, alpha, beta, m) -
                                        beta / (alpha + beta) *
                                        I1(R, am1, b, alpha, beta, m + 1))
    if b[i] > 0:
        bm1 = b[:i] + (b[i] - 1,) + b[i + 1:]
        result += b[i] / (2 * (alpha + beta)) * I1(R,
                                                   a, bm1, alpha, beta, m + 1)
    return result


def test_derivatives(R, a, b, alpha, beta, i):
    R = np.array(R)
    a = np.array(a)
    a[i] += 1
    dIdRi = 2 * alpha * I1(R, tuple(a), b, alpha, beta)
    a[i] -= 2
    if a[i] >= 0:
        dIdRi -= (a[i] + 1) * I1(R, tuple(a), b, alpha, beta)
    a[i] += 1
    dr = 0.001
    R[i] += 0.5 * dr
    dIdRi2 = I1(R, tuple(a), b, alpha, beta)
    R[i] -= dr
    dIdRi2 -= I1(R, tuple(a), b, alpha, beta)
    dIdRi2 /= -dr
    R[i] += 0.5 * dr
    return dIdRi, dIdRi2

class Gauss:
    """Normalised Gauss distribution

    from gauss import Gauss

    width=0.4
    gs = Gauss(width)

    for i in range(4):
        print 'Gauss(i)=',gs.get(i)
    """
    def __init__(self, width=0.08):
        self.dtype = float
        self.set_width(width)
        
    def get(self, x, x0=0):
        return self.norm * np.exp(-((x-x0) * self.wm1)**2)
    
    def set_width(self, width=0.08):
        self.norm = 1. / width / np.sqrt(2 * np.pi)
        self.wm1 = np.sqrt(.5) / width

class Lorentz:
    """Normalised Lorentz distribution"""
    def __init__(self, width=0.08):
        self.dtype = float
        self.set_width(width)
        
    def get(self, x, x0=0):
        return self.norm / ((x-x0)**2 + self.width2)
    
    def set_width(self, width=0.08):
        self.norm = width / np.pi
        self.width2 = width**2

