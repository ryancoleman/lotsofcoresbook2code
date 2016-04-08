# Copyright (C) 2006  CAMP
# Please see the accompanying LICENSE file for further information.

from math import pi, log, sqrt

import numpy as np

from gpaw.utilities import _fact

"""Fourier filtering

This module is an implementation of this Fourier filtering scheme:

*A general and efficient pseudopotential Fourier filtering scheme for
real space methods using mask functions*, Maxim Tafipolsky, Rochus
Schmid, J Chem Phys. 2006 May 7;124:174102.

Only difference is that we use a gaussian for the mask function.  The
filtering is used for the projector functions and for the zero
potential."""

# 3D-Fourier transform:
#
#                 /     _ _
# ~         ^    |  _  iq.r          ^
# f (q) Y  (q) = | dr e    f (r) Y  (r)
#  l     lm      |          l     lm
#               /
#
# Radial part:
#
#                  /
#   ~        __ l |  2
#   f (q) = 4||i  | r dr j (qr) f (r)
#    l            |       l      l
#                /
#

# XXX use fast bessel transform !!!

class Filter:
    """Mask-function Fourier filter"""
    
    def __init__(self, r_g, dr_g, gcut, h):
        """Construct filter.

        The radial grid is defined by r(g) and dr/dg(g) (`r_g` and
        `dr_g`), `gcut` is the cutoff grid point, and `h` is the target
        grid spacing used in the calculation."""

        self.gcut = gcut
        rcut = r_g[gcut]
        
        N = 200
        self.r_g = r_g = r_g[:gcut].copy()  # will be modified later!
        self.dr_g = dr_g[:gcut]

        # Matrices for Bessel transform:
        q1 = 5 * pi / h / N
        self.q_i = q_i = q1 * np.arange(N)
        self.c = sqrt(2 * q1 / pi) 
        self.sinqr_ig = np.sin(q_i[:, None] * r_g) * self.c
        self.cosqr_ig = np.cos(q_i[:, None] * r_g) * self.c

        # Cutoff function:
        qmax = pi / h
        alpha = 1.1
        qcut = qmax / alpha
        icut = 1 + int(qcut / q1)
        beta = 5 * log(10) / (alpha - 1.0)**2
        self.cut_i = np.ones(N)
        self.cut_i[icut:] = np.exp(
            -np.clip(beta * (q_i[icut:] / qcut - 1.0)**2, 0, 400))
        # self.cut_i[icut:] = np.exp(
        #     -np.clip(0, 400, beta * (q_i[icut:] / qcut - 1.0)**2))

        # Mask function:
        gamma = 3 * log(10) / rcut**2
        self.m_g = np.exp(-gamma * r_g**2)
        
        # We will need to divide by these two!  Remove zeros:
        q_i[0] = 1.0
        r_g[0] = 1.0

    def filter(self, f_g, l=0):
        """Filter radial function.

        The function to be filtered is::

          f(r)     ^
          ---- Y  (r)
           r    lm
           
        Output is::

                l     ^
          g(r) r  Y  (r),
                   lm

        where the filtered radial part ``g(r)`` is returned."""
        
        r_g = self.r_g
        q_i = self.q_i

        fdrim_g = f_g[:self.gcut] * self.dr_g / self.m_g / r_g

        #         sin(x)
        # j (x) = ------,
        #  0        x
        #
        #         sin(x)   cos(x)
        # j (x) = ------ - ------,
        #  1         2       x
        #           x
        #
        #           3    1            3
        # j (x) = (--- - -) sin(x) - --- cos(x),
        #  2         3   x             2
        #           x                 x
        #
        #          15    6             15   1 
        # j (x) = (-- - ---) sin(x) - (-- - -) cos(x).
        #  3        4     2             3   x 
        #          x     x             x      
        #


        if l == 0:
            fq_i = np.dot(self.sinqr_ig, fdrim_g * r_g) * self.cut_i
            fr_g = np.dot(fq_i, self.sinqr_ig)
        elif l == 1:
            fq_i = np.dot(self.sinqr_ig, fdrim_g) / q_i
            fq_i -= np.dot(self.cosqr_ig, r_g * fdrim_g)
            fq_i[0] = 0.0
            fq_i *= self.cut_i
            fr_g = np.dot(fq_i / q_i, self.sinqr_ig) / r_g
            fr_g -= np.dot(fq_i, self.cosqr_ig)
        elif l == 2:
            fq_i = 3 * np.dot(self.sinqr_ig, fdrim_g / r_g) / q_i**2
            fq_i -= np.dot(self.sinqr_ig, fdrim_g * r_g)
            fq_i -= 3 * np.dot(self.cosqr_ig, fdrim_g) / q_i
            fq_i[0] = 0.0
            fq_i *= self.cut_i
            fr_g = 3 * np.dot(fq_i / q_i**2, self.sinqr_ig) / r_g**2
            fr_g -= np.dot(fq_i, self.sinqr_ig)
            fr_g -= 3 * np.dot(fq_i / q_i, self.cosqr_ig) / r_g
        elif l == 3: ### This should be tested
            fq_i = 15 * np.dot(self.sinqr_ig, fdrim_g / r_g**2) / q_i**3
            fq_i -= 6 * np.dot(self.sinqr_ig, fdrim_g) / q_i
            fq_i -= 15 * np.dot(self.cosqr_ig, fdrim_g / r_g) / q_i**2
            fq_i += np.dot(self.cosqr_ig, r_g * fdrim_g)
            fq_i[0] = 0.0
            fq_i *= self.cut_i
            fr_g = 15 * np.dot(fq_i / q_i**3, self.sinqr_ig) / r_g**3
            fr_g -= 6 * np.dot(fq_i / q_i, self.sinqr_ig) / r_g
            fr_g -= 15 * np.dot(fq_i / q_i**2, self.cosqr_ig) / r_g**2
            fr_g += np.dot(fq_i, self.cosqr_ig)

        else:
            raise NotImplementedError
    
        a_g = np.zeros(len(f_g))
        a_g[:self.gcut] = fr_g * self.m_g / r_g**(l + 1)
        
        #            n 
        #           2 n!     n
        # j (x) = --------- x   for  x << 1.
        #  n      (2n + 1)! 
        #
        # This formula is used for finding the value of
        #
        #       -l
        # f(r) r    for r -> 0
        #
        c = 2.0**l * _fact[l] / _fact[2 * l + 1] * self.c
        a_g[0] = np.dot(fq_i, q_i**(l + 1)) * c

        return a_g

if __name__ == '__main__':
    rc = 1.1
    gamma = 1.95
    rc2 = rc * gamma
    M = 300
    beta = 0.3
    gcut = 1 + int(M * rc / (beta + rc))
    g_g = np.arange(M)
    r_g = beta * g_g / (M - g_g)
    drdg_g = beta * M / (M - g_g)**2

    x_g = r_g / rc
    p_g = 1 - x_g**2 * (3 - 2 * x_g)
    p_g[gcut:] = 0.0
    #p_g = np.exp(-np.clip(5.0 * r_g**2, 0, 400))

    h = 0.4
    f = Filter(r_g, drdg_g, rc2, h)
    pf0_g = f.filter(p_g)
    pf1_g = f.filter(p_g * r_g**1, 1)
    pf2_g = f.filter(p_g * r_g**2, 2)

    if 0:
        for i in range(200):
            print(5 * pi / h * i / 200, pf0_g[i], pf1_g[i], pf2_g[i])
    if 1:
        for r, p, pf0, pf1, pf2 in zip(r_g, p_g, pf0_g, pf1_g, pf2_g):
            print(r, p, pf0, pf1, pf2)
            if r > rc2:
                break
