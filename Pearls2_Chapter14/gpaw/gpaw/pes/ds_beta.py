from math import sin, cos, pi, sqrt

import numpy as np

from ase.units import Bohr, Hartree, alpha
from gpaw.fd_operators import Gradient
from gpaw.utilities.gl_quadrature import GaussLegendre
from gpaw.pes import ds_prefactor
from gpaw.pes.state import H1s
from gpaw.pes.continuum import PlaneWave

# debug
# from gpaw.mpi import rank


class CrossSectionBeta:

    def __init__(self,
                 initial=None,
                 final=None,
                 r0=[0, 0, 0],  # center of mass vector
                 form='L',
                 ngauss=8):

        self.initial = initial
        self.r0 = np.array(r0) / Bohr
        self.form = form
        if final is None:
            self.final = PlaneWave(initial.gd)
        else:
            self.final = final

        self.Ekin = None

        # Gauss-Legendre weights and abscissas
        self.gl = {}
        self.gl['x'] = GaussLegendre(-1., 1., ngauss)
        self.gl['phi'] = GaussLegendre(0, 2 * pi, ngauss)
        self.gl['psi'] = self.gl['phi']
        self.angle = {}

        # sin and cos of the magic angle (54.7 deg)
        self.costhm = 1. / sqrt(3)
        self.sinthm = sqrt(2. / 3.)

    def calculate(self, Ekin):
        """Calculate the necessary overlaps."""

        Ekin = Ekin / Hartree
        if self.Ekin == Ekin:
            return
        self.Ekin = Ekin

        # photoelectron momentum
        self.k = sqrt(2 * self.Ekin)

        for angle in ['x', 'phi', 'psi']:
            self.angle[angle] = self.gl[angle].get_x()[0]
        self.T20, self.T2m = self.gi_x()

        # we need the average
        self.T20 /= 8 * pi ** 2
        self.T2m /= 8 * pi ** 2

    def get_omega(self):
        """Return the necessary photon energy."""
        return self.Ekin - self.initial.get_energy() / Hartree

    def get_beta(self, Ekin=None):
        """Return the asymmetry parameter.

        E: photoelectron kinetic energy [eV]
        """
        if Ekin is not None:
            self.calculate(Ekin)
        return self.T20 / self.T2m - 1.

    def get_ds(self, Ekin=None, units='Mb'):
        """Return the total cross section.

        Ekin: photoelectron kinetic energy [eV]
        units:
        'Mb', 1 Mb = 1.e-22 m**2
        'Ang', 1 A**2 = 1.e-20 m**2
        'a.u.', 1 a_0**2 = 2.8e-21 m**2
        as output units
        """
        if Ekin is not None:
            self.calculate(Ekin)
        try:
            pre = ds_prefactor[units]
        except KeyError:
            print('Unknown units: >' + units + '<')
            raise

#        me_c =  self.initial.get_me_c(np.array([0., 0., self.k]), self.form)
#        T2mana = np.abs(np.dot(me_c,me_c)) / 3.
#        print "T2m:", T2mana, self.T2m

        omega = self.get_omega()

        # integration over momentum agles
        pre *= self.k * 4 * pi

#        print omega, self.initial.get_ds(self.k, omega, self.form), \
#            (self.k * 4 * pi * (2 * pi)**2 / 137.0359895 * self.T2m / omega)

        return pre * ((2 * pi) ** 2 * alpha * self.T2m / omega)

    def gauss_integrate(self, angle, function):
        T20 = 0.
        T2m = 0.
        gl = self.gl[angle]
        for x, w in zip(gl.get_x(), gl.get_w()):
#            print angle, x, w
            self.angle[angle] = x
            t20, t2m = function()
            T20 += w * t20
            T2m += w * t2m
#        print "<gauss_integrate> angle=", angle, T2m
        return T20, T2m

    def gi_x(self):
        """Gauss integrate over x=cos(theta)"""
        return self.gauss_integrate('x', self.gi_phi)

    def gi_phi(self):
        """Gauss integrate over phi"""
        return self.gauss_integrate('phi', self.gi_psi)

    def gi_psi(self):
        """Gauss integrate over psi"""
        return self.gauss_integrate('psi', self.integrand)

    def integrand(self):

        # polarisation in the direction of vk
        costh = self.angle['x']
        sinth = sqrt(1. - costh ** 2)
        sinphi = sin(self.angle['phi'])
        cosphi = cos(self.angle['phi'])
        eps0 = np.array([sinth * cosphi,
                         sinth * sinphi,
                         costh])
        vk = self.k * eps0

        # polarisation at the magic angle
        costhm = self.costhm
        sinthm = self.sinthm
        sinpsi = sin(self.angle['psi'])
        cospsi = cos(self.angle['psi'])
        epsm = np.array([sinthm * (cosphi * sinpsi * costh +
                                   sinphi * cospsi) +
                         costhm * cosphi * sinth,
                         sinthm * (sinphi * sinpsi * costh -
                                   cosphi * cospsi) +
                         costhm * sinphi * sinth,
                         costhm * costh - sinthm * sinth * sinpsi])

        # initial and final state on the grid
        initial_G = self.initial.get_grid()
        final_G = self.final.get_grid(vk, self.r0)
        ini_analyt = H1s(self.initial.gd, self.r0)

        gd = self.initial.gd
        if self.form == 'L':
            if_G = initial_G * final_G
            omega = self.get_omega()
            if 0:
                me_c = []
                for c in range(3):
                    xyz_G = ((np.arange(gd.n_c[c]) + gd.beg_c[c]) * gd.h_c[c]
                             - self.r0[c])
                    shape = [1, 1, 1]
                    shape[c] = -1
                    xyz_G.shape = tuple(shape)
                    np.resize(xyz_G, gd.n_c)
                    me_c.append(gd.integrate(if_G * xyz_G))
                me_c = np.array(me_c) * omega
            else:
                me_c = gd.calculate_dipole_moment(if_G)
                me_c += self.r0 * gd.integrate(if_G)
                me_c *= -omega
        elif self.form == 'V':
            dtype = final_G.dtype
            phase_cd = np.ones((3, 2), dtype)
            if not hasattr(gd, 'ddr'):
                gd.ddr = [Gradient(gd, c, dtype=dtype).apply for c in range(3)]
            dfinal_G = gd.empty(dtype=dtype)
            me_c = np.empty(3, dtype=dtype)
            for c in range(3):
                gd.ddr[c](final_G, dfinal_G, phase_cd)
                me_c[c] = gd.integrate(initial_G * dfinal_G)
        else:
            raise NotImplementedError

        if 0:
            omega = self.get_omega()
            me_analyt = ini_analyt.get_me_c(vk, self.form)[0].imag
            me = me_c[0].imag

            def ds(me):
                return self.k / omega * me ** 2

            print(omega, ds(me_analyt), ds(me), me_analyt, me)
#        print 'analyt', self.initial.get_me_c(vk, self.form)
#        print 'num', me_c
#        print 'analyt/num', self.initial.get_me_c(vk, self.form) / me_c

        # return the squared matrix elements
        T2 = []
        for eps in [eps0, epsm]:
            me = np.dot(eps, me_c)
#            print "eps, T2:", eps, (me * me.conj()).real
            T2.append((me * me.conj()).real)
#        print vk, T2
        return T2[0], T2[1]
