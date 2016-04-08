from __future__ import division

from math import pi

import numpy as np

from gpaw.spline import Spline
from gpaw.utilities import hartree, divrl, _fact as fac


def radial_grid_descriptor(eq, **kwargs):
    if eq == 'r=d*i':
        assert int(kwargs['istart']) == 0
        return EquidistantRadialGridDescriptor(float(kwargs['d']),
                                               int(kwargs['iend']) + 1)
    if eq == 'r=a*i/(n-i)':
        beta = float(kwargs['a'])
        ng = int(kwargs['n'])
        return AERadialGridDescriptor(beta / ng, 1.0 / ng, ng)
    if eq == 'r=a*i/(1-b*i)':
        a = float(kwargs['a'])
        b = float(kwargs['b'])
        N = int(kwargs['n'])
        return AERadialGridDescriptor(a, b, N)
    if eq == 'r=a*(exp(d*i)-1)':
        a = float(kwargs['a'])
        d = float(kwargs['d'])
        N = int(kwargs['iend']) + 1
        return AbinitRadialGridDescriptor(a, d, N)
    raise ValueError('Unknown grid: ' + eq)


def fsbt(l, f_g, r_g, G_k):
    """Fast spherical Bessel transform.

    Returns::
       
          oo
         / 2
         |r dr j (Gr) f(r),
         /      l
          0

    using l+1 fft's."""

    N = (len(G_k) - 1) * 2
    f_k = 0.0
    F_g = f_g * r_g
    for n in range(l + 1):
        f_k += (r_g[1] * (1j)**(l + 1 - n) *
                fac[l + n] / fac[l - n] / fac[n] / 2**n *
                np.fft.rfft(F_g, N)).real * G_k**(l - n)
        F_g[1:] /= r_g[1:]

    f_k[1:] /= G_k[1:]**(l + 1)
    if l == 0:
        f_k[0] = np.dot(r_g, f_g * r_g) * r_g[1]
    return f_k


class RadialGridDescriptor:
    ndim = 1  # dimension of ndarrays

    def __init__(self, r_g, dr_g, default_spline_points=25):
        """Grid descriptor for radial grid."""
        self.r_g = r_g
        self.dr_g = dr_g
        self.N = len(r_g)
        self.dv_g = 4 * pi * r_g**2 * dr_g
        self.default_spline_points = default_spline_points

    def __len__(self):
        return self.N

    def zeros(self, x=()):
        a_xg = self.empty(x)
        a_xg[:] = 0
        return a_xg

    def empty(self, x=()):
        if isinstance(x, int):
            x = (x,)
        return np.zeros(x + (self.N,))

    def integrate(self, a_xg, n=0):
        assert n >= -2
        return np.dot(a_xg[..., 1:],
                      (self.r_g**(2 + n) * self.dr_g)[1:]) * (4 * pi)

    def integrate_yukawa(self, n1, n2, l, gamma):
        """Integrate two densities n1 and n2 with yukawa interaction."""
        from scipy.special import iv, kv
        r = self.r_g
        dr = self.dr_g
        k_rgamma = kv(l + 0.5, r * gamma)      # K(>)
        i_rgamma = iv(l + 0.5, r * gamma)      # I(<)
        k_rgamma[0] = kv(l + 0.5, r[1] * gamma * 1e-5)
        matrix_ik = np.outer(n1 * dr, n2 * dr)
        len_vec = len(k_rgamma)
        for i in xrange(len_vec):
            k_rgi = k_rgamma[i]
            for k in xrange(i):
                modified_bessels = i_rgamma[k] * k_rgi
                matrix_ik[i, k] *= modified_bessels
                matrix_ik[k, i] *= modified_bessels
            matrix_ik[i, i] *= i_rgamma[i] * k_rgi
        return matrix_ik.sum()

    def derivative(self, n_g, dndr_g=None):
        """Finite-difference derivative of radial function."""
        if dndr_g is None:
            dndr_g = self.empty()
        dndr_g[0] = n_g[1] - n_g[0]
        dndr_g[1:-1] = 0.5 * (n_g[2:] - n_g[:-2])
        dndr_g[-1] = n_g[-1] - n_g[-2]
        dndr_g /= self.dr_g
        return dndr_g

    def derivative2(self, a_g, b_g):
        """Finite-difference derivative of radial function.

        For an infinitely dense grid, this method would be identical
        to the `derivative` method."""
        
        c_g = a_g / self.dr_g
        b_g[0] = 0.5 * c_g[1] + c_g[0]
        b_g[1] = 0.5 * c_g[2] - c_g[0]
        b_g[1:-1] = 0.5 * (c_g[2:] - c_g[:-2])
        b_g[-2] = c_g[-1] - 0.5 * c_g[-3]
        b_g[-1] = -c_g[-1] - 0.5 * c_g[-2]

    def interpolate(self, f_g, r_x):
        from scipy.interpolate import InterpolatedUnivariateSpline
        return InterpolatedUnivariateSpline(self.r_g, f_g)(r_x)
        
    def fft(self, fr_g, l=0, N=None):
        """Fourier transform.

        Returns G and f(G) arrays::
           
                                          _ _
               l    ^    / _         ^   iG.r
          f(G)i Y  (G) = |dr f(r)Y  (r) e    .
                 lm      /        lm
        """

        if N is None:
            N = 2**13

        assert N % 2 == 0

        r_x = np.linspace(0, self.r_g[-1], N)
        f_x = self.interpolate(fr_g, r_x)
        f_x[1:] /= r_x[1:]
        f_x[0] = f_x[1]
        G_k = np.linspace(0, pi / r_x[1], N // 2 + 1)
        f_k = 4 * pi * fsbt(l, f_x, r_x, G_k)
        return G_k, f_k

    def filter(self, f_g, rcut, Gcut, l=0):
        Rcut = 100.0
        N = 1024 * 8
        r_x = np.linspace(0, Rcut, N, endpoint=False)
        h = Rcut / N

        alpha = 4.0 / rcut**2
        mcut = np.exp(-alpha * rcut**2)
        r2_x = r_x**2
        m_x = np.exp(-alpha * r2_x)
        for n in range(2):
            m_x -= (alpha * (rcut**2 - r2_x))**n * (mcut / fac[n])
        xcut = int(np.ceil(rcut / r_x[1]))
        m_x[xcut:] = 0.0

        G_k = np.linspace(0, pi / h, N // 2 + 1)

        from scipy.interpolate import InterpolatedUnivariateSpline
        if l < 2:
            f_x = InterpolatedUnivariateSpline(self.r_g, f_g)(r_x)
        else:
            a_g = f_g.copy()
            a_g[1:] /= self.r_g[1:]**(l - 1)
            f_x = InterpolatedUnivariateSpline(
                self.r_g, a_g)(r_x) * r_x**(l - 1)

        f_x[:xcut] /= m_x[:xcut]
        f_k = fsbt(l, f_x, r_x, G_k)
        kcut = int(Gcut / G_k[1])
        f_k[kcut:] = 0.0
        ff_x = fsbt(l, f_k, G_k, r_x[:N // 2 + 1]) / pi * 2
        ff_x *= m_x[:N // 2 + 1]

        if l < 2:
            f_g = InterpolatedUnivariateSpline(
                r_x[:xcut + 1], ff_x[:xcut + 1])(self.r_g)
        else:
            ff_x[1:xcut + 1] /= r_x[1:xcut + 1]**(l - 1)
            f_g = InterpolatedUnivariateSpline(
                r_x[:xcut + 1], ff_x[:xcut + 1])(self.r_g) * self.r_g**(l - 1)
        f_g[self.ceil(rcut):] = 0.0

        return f_g

    def poisson(self, n_g, l=0):
        vr_g = self.zeros()
        nrdr_g = n_g * self.r_g * self.dr_g
        hartree(l, nrdr_g, self.r_g, vr_g)
        return vr_g

    def pseudize(self, a_g, gc, l=0, points=4):
        """Construct smooth continuation of a_g for g<gc.
        
        Returns (b_g, c_p[P-1]) such that b_g=a_g for g >= gc and::
        
                P-1      2(P-1-p)+l
            b = Sum c_p r
             g  p=0      g

        for g < gc+P.
        """
        assert isinstance(gc, int) and gc > 10
        
        r_g = self.r_g
        i = range(gc, gc + points)
        r_i = r_g[i]
        c_p = np.polyfit(r_i**2, a_g[i] / r_i**l, points - 1)
        b_g = a_g.copy()
        b_g[:gc] = np.polyval(c_p, r_g[:gc]**2) * r_g[:gc]**l
        return b_g, c_p[-1]

    def pseudize_normalized(self, a_g, gc, l=0, points=3):
        """Construct normalized smooth continuation of a_g for g<gc.
        
        Same as pseudize() with also this constraint::
        
            /        2  /        2
            | dr b(r) = | dr a(r)
            /           /
        """

        b_g = self.pseudize(a_g, gc, l, points)[0]
        c_x = np.empty(points + 1)
        gc0 = gc // 2
        x0 = b_g[gc0]
        r_g = self.r_g
        i = [gc0] + range(gc, gc + points)
        r_i = r_g[i]
        norm = self.integrate(a_g**2)

        def f(x):
            b_g[gc0] = x
            c_x[:] = np.polyfit(r_i**2, b_g[i] / r_i**l, points)
            b_g[:gc] = np.polyval(c_x, r_g[:gc]**2) * r_g[:gc]**l
            return self.integrate(b_g**2) - norm

        from scipy.optimize import fsolve
        fsolve(f, x0)
        return b_g, c_x[-1]

    def jpseudize(self, a_g, gc, l=0, points=4):
        """Construct spherical Bessel function continuation of a_g for g<gc.
        
        Returns (b_g, b(0)/r^l) such that b_g=a_g for g >= gc and::
        
                P-2
            b = Sum c_p j (q r )
             g  p=0      l  p g

        for g < gc+P, where.
        """

        from scipy.special import sph_jn
        from scipy.optimize import brentq

        if a_g[gc] == 0:
            return self.zeros(), 0.0

        assert isinstance(gc, int) and gc > 10

        zeros_l = [[1, 2, 3, 4, 5, 6],
                   [1.430, 2.459, 3.471, 4.477, 5.482, 6.484],
                   [1.835, 2.895, 3.923, 4.938, 5.949, 6.956],
                   [2.224, 3.316, 4.360, 5.387, 6.405, 7.418]]
        
        # Logarithmic derivative:
        ld = np.dot([-1 / 60, 3 / 20, -3 / 4, 0, 3 / 4, -3 / 20, 1 / 60],
                    a_g[gc - 3:gc + 4]) / a_g[gc] / self.dr_g[gc]
        
        rc = self.r_g[gc]

        def f(q):
            j, dj = (y[-1] for y in sph_jn(l, q * rc))
            return dj * q - j * ld

        j_pg = self.empty(points - 1)
        q_p = np.empty(points - 1)

        zeros = zeros_l[l]
        if rc * ld > l:
            z1 = zeros[0]
            zeros = zeros[1:]
        else:
            z1 = 0
        x = 0.01
        for p, z2 in enumerate(zeros[:points - 1]):
            q = brentq(f, z1 * pi / rc + x, z2 * pi / rc - x)
            j_pg[p] = [sph_jn(l, q * r)[0][-1] for r in self.r_g]
            q_p[p] = q
            z1 = z2

        C_dg = [[0, 0, 0, 1, 0, 0, 0],
                [1 / 90, -3 / 20, 3 / 2, -49 / 18, 3 / 2, -3 / 20, 1 / 90],
                [1 / 8, -1, 13 / 8, 0, -13 / 8, 1, -1 / 8],
                [-1 / 6, 2, -13 / 2, 28 / 3, -13 / 2, 2, -1 / 6]][:points - 1]
        c_p = np.linalg.solve(np.dot(C_dg, j_pg[:, gc - 3:gc + 4].T),
                              np.dot(C_dg, a_g[gc - 3:gc + 4]))
        b_g = a_g.copy()
        b_g[:gc + 2] = np.dot(c_p, j_pg[:, :gc + 2])
        return b_g, np.dot(c_p, q_p**l) * 2**l * fac[l] / fac[2 * l + 1]

    def plot(self, a_g, n=0, rc=4.0, show=False):
        import matplotlib.pyplot as plt
        r_g = self.r_g[:len(a_g)]
        if n < 0:
            r_g = r_g[1:]
            a_g = a_g[1:]
        plt.plot(r_g, a_g * r_g**n)
        plt.axis(xmax=rc)
        if show:
            plt.show()

    def floor(self, r):
        return np.floor(self.r2g(r)).astype(int)
    
    def round(self, r):
        return np.around(self.r2g(r)).astype(int)
    
    def ceil(self, r):
        return np.ceil(self.r2g(r)).astype(int)

    def spline(self, a_g, rcut=None, l=0, points=None):
        if points is None:
            points = self.default_spline_points

        if rcut is None:
            g = len(a_g) - 1
            while a_g[g] == 0.0:
                g -= 1
            rcut = self.r_g[g + 1]

        b_g = a_g.copy()
        N = len(b_g)
        if l > 0:
            b_g = divrl(b_g, l, self.r_g[:N])
            
        r_i = np.linspace(0, rcut, points + 1)
        g_i = np.clip((self.r2g(r_i) + 0.5).astype(int), 1, N - 2)

        r1_i = self.r_g[g_i - 1]
        r2_i = self.r_g[g_i]
        r3_i = self.r_g[g_i + 1]
        x1_i = (r_i - r2_i) * (r_i - r3_i) / (r1_i - r2_i) / (r1_i - r3_i)
        x2_i = (r_i - r1_i) * (r_i - r3_i) / (r2_i - r1_i) / (r2_i - r3_i)
        x3_i = (r_i - r1_i) * (r_i - r2_i) / (r3_i - r1_i) / (r3_i - r2_i)
        b1_i = b_g[g_i - 1]
        b2_i = b_g[g_i]
        b3_i = b_g[g_i + 1]
        b_i = b1_i * x1_i + b2_i * x2_i + b3_i * x3_i
        return Spline(l, rcut, b_i)


class EquidistantRadialGridDescriptor(RadialGridDescriptor):
    def __init__(self, h, N=1000, h0=0.0):
        """Equidistant radial grid descriptor.

        The radial grid is r(g) = h0 + g*h,  g = 0, 1, ..., N - 1."""

        RadialGridDescriptor.__init__(self, h * np.arange(N) + h0, h)

    def r2g(self, r):
        return (r - self.r_g[0]) / (self.r_g[1] - self.r_g[0])

    def spline(self, a_g, rcut=None, l=0):
        assert rcut is None
        b_g = a_g.copy()
        if l > 0:
            b_g = divrl(b_g, l, self.r_g[:len(a_g)])
        return Spline(l, self.r_g[len(a_g) - 1], b_g)

        
class AERadialGridDescriptor(RadialGridDescriptor):
    def __init__(self, a, b, N=1000, default_spline_points=25):
        """Radial grid descriptor for all-electron calculation.

        The radial grid is::

                     a g
            r(g) = -------,  g = 0, 1, ..., N - 1
                   1 - b g
        """

        self.a = a
        self.b = b
        g = np.arange(N)
        r_g = self.a * g / (1 - self.b * g)
        dr_g = (self.b * r_g + self.a)**2 / self.a
        RadialGridDescriptor.__init__(self, r_g, dr_g, default_spline_points)
                                      
    def r2g(self, r):
        # return r / (r * self.b + self.a)
        # Hack to preserve backwards compatibility:
        ng = 1.0 / self.b
        beta = self.a / self.b
        return ng * r / (beta + r)

    def xml(self, id='grid1'):
        if abs(self.N - 1 / self.b) < 1e-5:
            return (('<radial_grid eq="r=a*i/(n-i)" a="%r" n="%d" ' +
                     'istart="0" iend="%d" id="%s"/>\n') %
                    (self.a * self.N, self.N, self.N - 1, id))
        return (('<radial_grid eq="r=a*i/(1-b*i)" a="%r" b="%r" n="%d" ' +
                 'istart="0" iend="%d" id="%s"/>\n') %
                (self.a, self.b, self.N, self.N - 1, id))

    def d2gdr2(self):
        return -2 * self.a * self.b / (self.b * self.r_g + self.a)**3

        
class AbinitRadialGridDescriptor(RadialGridDescriptor):
    def __init__(self, a, d, N=1000, default_spline_points=25):
        """Radial grid descriptor for Abinit calculations.

        The radial grid is::

                      dg
            r(g) = a(e   - 1),  g = 0, 1, ..., N - 1
        """

        self.a = a
        self.d = d
        g = np.arange(N)
        r_g = a * (np.exp(d * g) - 1)
        dr_g = (r_g + a) * d
        RadialGridDescriptor.__init__(self, r_g, dr_g, default_spline_points)
                                      
    def r2g(self, r):
        return np.log(r / self.a + 1) / self.d
