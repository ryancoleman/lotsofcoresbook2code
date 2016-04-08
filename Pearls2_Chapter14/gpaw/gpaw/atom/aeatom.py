# -*- coding: utf-8 -*-
from __future__ import print_function
import copy
import sys
from math import pi

import numpy as np
from numpy.linalg import eigh
from scipy.special import gamma
# from scipy.linalg import solve_banded
import ase.units as units
from ase.data import atomic_numbers, atomic_names, chemical_symbols

from gpaw.xc import XC
from gpaw.gaunt import make_gaunt
from gpaw.atom.configurations import configurations
from gpaw.atom.radialgd import AERadialGridDescriptor


# Velocity of light in atomic units:
c = 2 * units._hplanck / (units._mu0 * units._c * units._e**2)

# Colors for s, p, d, f, g:
colors = 'krgbycmmmmm'


class GaussianBasis:
    def __init__(self, l, alpha_B, rgd, eps=1.0e-7):
        """Guassian basis set for spherically symmetric atom.

        l: int
            Angular momentum quantum number.
        alpha_B: ndarray
            Exponents.
        rgd: GridDescriptor
            Grid descriptor.
        eps: float
            Cutoff for eigenvalues of overlap matrix."""
        
        self.l = l
        self.alpha_B = alpha_B
        self.rgd = rgd
        self.eps = eps

        A_BB = np.add.outer(alpha_B, alpha_B)
        M_BB = np.multiply.outer(alpha_B, alpha_B)

        # Overlap matrix:
        S_BB = (2 * M_BB**0.5 / A_BB)**(l + 1.5)

        # Kinetic energy matrix:
        T_BB = 2**(l + 2.5) * M_BB**(0.5 * l + 0.75) / gamma(l + 1.5) * (
            gamma(l + 2.5) * M_BB / A_BB**(l + 2.5) -
            0.5 * (l + 1) * gamma(l + 1.5) / A_BB**(l + 0.5) +
            0.25 * (l + 1) * (2 * l + 1) * gamma(l + 0.5) / A_BB**(l + 0.5))

        # Derivative matrix:
        D_BB = 2**(l + 2.5) * M_BB**(0.5 * l + 0.75) / gamma(l + 1.5) * (
            0.5 * (l + 1) * gamma(l + 1) / A_BB**(l + 1) -
            gamma(l + 2) * alpha_B / A_BB**(l + 2))

        # 1/r matrix:
        K_BB = 2**(l + 2.5) * M_BB**(0.5 * l + 0.75) / gamma(l + 1.5) * (
            0.5 * gamma(l + 1) / A_BB**(l + 1))

        # Find set of linearly independent functions.
        # We have len(alpha_B) gaussians (index B) and self.nbasis
        # linearly independent functions (index b).
        s_B, U_BB = eigh(S_BB)
        self.nbasis = int((s_B > eps).sum())

        Q_Bb = np.dot(U_BB[:, -self.nbasis:],
                      np.diag(s_B[-self.nbasis:]**-0.5))
        
        self.T_bb = np.dot(np.dot(Q_Bb.T, T_BB), Q_Bb)
        self.D_bb = np.dot(np.dot(Q_Bb.T, D_BB), Q_Bb)
        self.K_bb = np.dot(np.dot(Q_Bb.T, K_BB), Q_Bb)

        r_g = rgd.r_g
        # Avoid errors in debug mode from division by zero:
        old_settings = np.seterr(divide='ignore')
        self.basis_bg = (np.dot(
                Q_Bb.T,
                (2 * (2 * alpha_B[:, None])**(l + 1.5) /
                 gamma(l + 1.5))**0.5 *
                np.exp(-np.multiply.outer(alpha_B, r_g**2))) * r_g**l)
        np.seterr(**old_settings)
        
    def __len__(self):
        return self.nbasis

    def expand(self, C_xb):
        return np.dot(C_xb, self.basis_bg)

    def calculate_potential_matrix(self, vr_g):
        vr2dr_g = vr_g * self.rgd.r_g * self.rgd.dr_g
        V_bb = np.inner(self.basis_bg[:, 1:],
                        self.basis_bg[:, 1:] * vr2dr_g[1:])
        return V_bb


def coefs(rgd, l, vr_g, e, scalar_relativistic=False, Z=None):
    r_g = rgd.r_g

    x0_g = 2 * (e * r_g - vr_g) * r_g
    x1_g = 2 * (l + 1) * r_g / rgd.dr_g + r_g**2 * rgd.d2gdr2()
    x2_g = r_g**2 / rgd.dr_g**2
    x = 1.0
    
    if scalar_relativistic:
        x = (1 + l + l**2 - (Z / c)**2)**0.5 - l
        r_g = r_g.copy()
        r_g[0] = 1.0
        v_g = vr_g / r_g
        M_g = 1 + (e - v_g) / (2 * c**2)
        kappa_g = (rgd.derivative(vr_g) - v_g) / r_g / (2 * c**2 * M_g)
        x0_g = (2 * M_g * (e * r_g - vr_g) * r_g +
                (l + x - 1) * kappa_g * r_g +
                (l + x) * (l + x - 1) - l * (l + 1))
        x1_g = (2 * (l + x) * r_g / rgd.dr_g +
                r_g**2 * rgd.d2gdr2() +
                r_g**2 * kappa_g / rgd.dr_g)

    cm1_g = x2_g - x1_g / 2
    c0_g = x0_g - 2 * x2_g
    cp1_g = x2_g + x1_g / 2
    
    return cm1_g, c0_g, cp1_g, x
    

class Channel:
    def __init__(self, l, s=0, f_n=(), basis=None):
        self.l = l
        self.s = s
        self.basis = basis

        self.C_nb = None                       # eigenvectors
        self.e_n = None                        # eigenvalues
        self.f_n = np.array(f_n, dtype=float)  # occupation numbers
        self.phi_ng = None                     # wave functions
        
        self.name = 'spdfg'[l]
        self.solve2ok = False
        
    def solve(self, vr_g):
        """Diagonalize Schrödinger equation in basis set."""
        H_bb = self.basis.calculate_potential_matrix(vr_g)
        H_bb += self.basis.T_bb
        self.e_n, C_bn = eigh(H_bb)
        self.C_nb = C_bn.T
        self.phi_ng = self.basis.expand(self.C_nb[:len(self.f_n)])

    def solve2(self, vr_g, scalar_relativistic=False, Z=None):
        rgd = self.basis.rgd
        r_g = rgd.r_g
        l = self.l
        u_g = rgd.empty()
        self.solve2ok = True
        for n in range(len(self.f_n)):
            e = self.e_n[n]

            # Find classical turning point:
            x = vr_g * r_g + 0.5 * l * (l + 1) - e * r_g**2
            g0 = rgd.round(4.0)
            while x[g0] > 0:
                g0 -= 1

            iter = 0
            ok = False
            while True:
                du1dr, a = self.integrate_outwards(u_g, rgd, vr_g, g0, e,
                                                   scalar_relativistic, Z)
                u1 = u_g[g0]
                du2dr = self.integrate_inwards(u_g, rgd, vr_g, g0, e,
                                               scalar_relativistic, Z)
                u2 = u_g[g0]
                A = du1dr / u1 - du2dr / u2
                u_g[g0:] *= u1 / u2
                norm = rgd.integrate(u_g**2, -2) / (4 * pi)
                u_g /= norm**0.5
                a /= norm**0.5

                nodes = (u_g[:-1] * u_g[1:] < 0).sum()

                if abs(A) < 1e-5 and nodes == n:
                    ok = True
                    break

                if nodes > n:
                    e *= 1.2
                elif nodes < n:
                    e *= 0.8
                else:
                    e += 0.5 * A * u_g[g0]**2
                    if e > 0:
                        break
                
                iter += 1
                assert iter < 400, (n, l, e)
                
            if ok:
                self.e_n[n] = e
                self.phi_ng[n, 1:] = u_g[1:] / r_g[1:]
                self.phi_ng[n, 0] = a
            else:
                self.solve2ok = False
            
    def calculate_density(self, n=None):
        """Calculate density."""
        if n is None:
            n_g = 0.0
            for n, f in enumerate(self.f_n):
                n_g += f * self.calculate_density(n)
        else:
            n_g = self.phi_ng[n]**2 / (4 * pi)
        return n_g

    def calculate_kinetic_energy_density(self, n):
        """Calculate kinetic energy density."""
        phi_g = self.phi_ng[n]
        rgd = self.basis.rgd
        tau_g = rgd.derivative(phi_g)**2 / (8 * pi)
        if self.l > 0:
            tau_g[1:] += (self.l * (self.l + 1) *
                          (phi_g[1:] / rgd.r_g[1:])**2 / (8 * pi))
        return tau_g

    def get_eigenvalue_sum(self):
        f_n = self.f_n
        return np.dot(f_n, self.e_n[:len(f_n)])

    def integrate_outwards(self, u_g, rgd, vr_g, g0, e,
                           scalar_relativistic=False, Z=None, pt_g=None):
        l = self.l
        r_g = rgd.r_g

        cm1_g, c0_g, cp1_g, x = coefs(rgd, l, vr_g, e, scalar_relativistic, Z)

        # c_xg = np.zeros((3, g0 + 2))
        # c_xg[0, :2] = 1.0
        # c_xg[0, 2:] = cp1_g[1:g0 + 1]
        # c_xg[1, 1:-1] = c0_g[1:g0 + 1]
        # c_xg[2, :-2] = cm1_g[1:g0 + 1]

        # b_g = np.zeros(g0 + 2)
        if pt_g is not None:
            2 / 0  # Need to fix this:
            # b_g[2:] = -2 * pt_g[1:g0 + 1] * r_g[1:g0 + 1]**(1 - l)
            # a0 = pt_g[1] / r_g[1]**l / (vr_g[1] / r_g[1] - e)

        # b_g[:2] = [a0, a1]
        # a_g = solve_banded((2, 0), c_xg, b_g,
        #                    overwrite_ab=True, overwrite_b=True)
        
        g = 1
        agm1 = 1.0
        u_g[0] = 0.0
        ag = 1.0
        while True:
            u_g[g] = ag * r_g[g]**(l + x)
            agp1 = -(agm1 * cm1_g[g] + ag * c0_g[g]) / cp1_g[g]
            if g == g0:
                break
            g += 1
            agm1 = ag
            ag = agp1

        r = r_g[g0]
        dr = rgd.dr_g[g0]
        da = 0.5 * (agp1 - agm1)
        dudr = (l + x) * r**(l + x - 1) * ag + r**(l + x) * da / dr

        if l - 1 + x < 0:
            phi0 = (r_g[1] * 0.1)**(l - 1 + x)
        else:
            phi0 = 0.0**(l - 1 + x)
            
        return dudr, phi0

    def integrate_inwards(self, u_g, rgd, vr_g, g0, e,
                          scalar_relativistic=False, Z=None, gmax=None):
        l = self.l
        r_g = rgd.r_g

        cm1_g, c0_g, cp1_g, x = coefs(rgd, l, vr_g, e, scalar_relativistic, Z)

        cm1_g[:g0] = 1.0  # prevent division by zero
        c0_g /= -cm1_g
        cp1_g /= -cm1_g

        if gmax is None:
            gmax = len(u_g)

        g = gmax - 2
        agp1 = 1.0
        u_g[gmax - 1] = agp1 * r_g[gmax - 1]**(l + x)
        ag = np.exp(-(-2 * e)**0.5 * (r_g[gmax - 2] - r_g[gmax - 1]))

        while True:
            u_g[g] = ag * r_g[g]**(l + x)
            if ag > 1e50:
                u_g[g:] /= 1e50
                ag = ag / 1e50
                agp1 = agp1 / 1e50
            agm1 = agp1 * cp1_g[g] + ag * c0_g[g]
            if g == g0:
                break
            g -= 1
            agp1 = ag
            ag = agm1

        r = r_g[g]
        dr = rgd.dr_g[g]
        da = 0.5 * (agp1 - agm1)
        dudr = (l + x) * r**(l + x - 1) * ag + r**(l + x) * da / dr

        return dudr


class DiracChannel(Channel):
    def __init__(self, k, f_n, basis):
        l = (abs(2 * k + 1) - 1) // 2
        Channel.__init__(self, l, 0, f_n, basis)
        self.k = k
        self.j = abs(k) - 0.5
        self.c_nb = None  # eigenvectors (small component)

        self.name += '(%d/2)' % (2 * self.j)

    def solve(self, vr_g):
        """Solve Dirac equation in basis set."""
        nb = len(self.basis)
        V_bb = self.basis.calculate_potential_matrix(vr_g)
        H_bb = np.zeros((2 * nb, 2 * nb))
        H_bb[:nb, :nb] = V_bb
        H_bb[nb:, nb:] = V_bb - 2 * c**2 * np.eye(nb)
        H_bb[nb:, :nb] = -c * (-self.basis.D_bb.T + self.k * self.basis.K_bb)
        e_n, C_bn = eigh(H_bb)
        if self.k < 0:
            n0 = nb
        else:
            n0 = nb + 1
        self.e_n = e_n[n0:].copy()
        self.C_nb = C_bn[:nb, n0:].T.copy()  # large component
        self.c_nb = C_bn[nb:, n0:].T.copy()  # small component

    def calculate_density(self, n=None):
        """Calculate density."""
        if n is None:
            n_g = Channel.calculate_density(self)
        else:
            n_g = (self.basis.expand(self.C_nb[n])**2 +
                   self.basis.expand(self.c_nb[n])**2) / (4 * pi)
            if self.basis.l < 0:
                n_g[0] = n_g[1]
        return n_g

        
class AllElectronAtom:
    def __init__(self, symbol, xc='LDA', spinpol=False, dirac=False,
                 configuration=None,
                 log=None):
        """All-electron calculation for spherically symmetric atom.

        symbol: str (or int)
            Chemical symbol (or atomic number).
        xc: str
            Name of XC-functional.
        spinpol: bool
            If true, do spin-polarized calculation.  Default is spin-paired.
        dirac: bool
            Solve Dirac equation instead of Schrödinger equation.
        configuration: list
            Electronic configuration for symbol, format as in
            gpaw.atom.configurations
        log: stream
            Text output."""

        if isinstance(symbol, int):
            symbol = chemical_symbols[symbol]
        self.symbol = symbol
        self.Z = atomic_numbers[symbol]

        self.nspins = 1 + int(bool(spinpol))

        self.dirac = bool(dirac)

        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)
        else:
            self.configuration = None
        
        self.scalar_relativistic = False

        if isinstance(xc, str):
            self.xc = XC(xc)
        else:
            self.xc = xc

        self.fd = log or sys.stdout

        self.vr_sg = None  # potential * r
        self.n_sg = 0.0  # density
        self.rgd = None  # radial grid descriptor

        # Energies:
        self.ekin = None
        self.eeig = None
        self.eH = None
        self.eZ = None

        self.channels = None

        self.initialize_configuration(self.configuration)

        self.log('Z:              ', self.Z)
        self.log('Name:           ', atomic_names[self.Z])
        self.log('Symbol:         ', symbol)
        self.log('XC-functional:  ', self.xc.name)
        self.log('Equation:       ', ['Schrödinger', 'Dirac'][self.dirac])

        self.method = 'Gaussian basis-set'

    def log(self, *args, **kwargs):
        print(file=self.fd, *args, **kwargs)

    def initialize_configuration(self, configuration=None):
        self.f_lsn = {}

        if configuration is None:
            configuration = configurations[self.symbol][1]

        for n, l, f, e in configuration:
            
            if l not in self.f_lsn:
                self.f_lsn[l] = [[] for s in range(self.nspins)]
            if self.nspins == 1:
                self.f_lsn[l][0].append(f)
            else:
                # Use Hund's rule:
                f0 = min(f, 2 * l + 1)
                self.f_lsn[l][0].append(f0)
                self.f_lsn[l][1].append(f - f0)
                
        if 0:
            n = 2 + len(self.f_lsn[2][0])
            if self.f_lsn[0][0][n] == 2:
                self.f_lsn[0][0][n] = 1
                self.f_lsn[2][0][n - 3] += 1

    def add(self, n, l, df=+1, s=None):
        """Add (remove) electrons."""
        if s is None:
            if self.nspins == 1:
                s = 0
            else:
                self.add(n, l, 0.5 * df, 0)
                self.add(n, l, 0.5 * df, 1)
                return
            
        if l not in self.f_lsn:
            self.f_lsn[l] = [[] for x in range(self.nspins)]
            
        f_n = self.f_lsn[l][s]
        if len(f_n) < n - l:
            f_n.extend([0] * (n - l - len(f_n)))
        f_n[n - l - 1] += df

    def initialize(self, ngpts=2000, rcut=50.0,
                   alpha1=0.01, alpha2=None, ngauss=50,
                   eps=1.0e-7):
        """Initialize basis sets and radial grid.

        ngpts: int
            Number of grid points for radial grid.
        rcut: float
            Cutoff for radial grid.
        alpha1: float
            Smallest exponent for gaussian.
        alpha2: float
            Largest exponent for gaussian.
        ngauss: int
            Number of gaussians.
        eps: float
            Cutoff for eigenvalues of overlap matrix."""

        if alpha2 is None:
            alpha2 = 50.0 * self.Z**2

        # Use grid with r(0)=0, r(1)=a and r(ngpts)=rcut:
        a = 1 / alpha2**0.5 / 20
        b = (rcut - a * ngpts) / (rcut * ngpts)
        b = 1 / round(1 / b)
        self.rgd = AERadialGridDescriptor(a, b, ngpts)
        
        self.log('Grid points:     %d (%.5f, %.5f, %.5f, ..., %.3f, %.3f)' %
                 ((self.rgd.N,) + tuple(self.rgd.r_g[[0, 1, 2, -2, -1]])))

        # Distribute exponents between alpha1 and alpha2:
        alpha_B = alpha1 * (alpha2 / alpha1)**np.linspace(0, 1, ngauss)
        self.log('Exponents:       %d (%.3f, %.3f, ..., %.3f, %.3f)' %
                 ((ngauss,) + tuple(alpha_B[[0, 1, -2, -1]])))

        # Maximum l value:
        lmax = max(self.f_lsn.keys())

        self.channels = []
        nb_l = []
        if not self.dirac:
            for l in range(lmax + 1):
                basis = GaussianBasis(l, alpha_B, self.rgd, eps)
                nb_l.append(len(basis))
                for s in range(self.nspins):
                    self.channels.append(Channel(l, s, self.f_lsn[l][s],
                                                 basis))
        else:
            for K in range(1, lmax + 2):
                leff = (K**2 - (self.Z / c)**2)**0.5 - 1
                basis = GaussianBasis(leff, alpha_B, self.rgd, eps)
                nb_l.append(len(basis))
                for k, l in [(-K, K - 1), (K, K)]:
                    if l > lmax:
                        continue
                    f_n = self.f_lsn[l][0]
                    j = abs(k) - 0.5
                    f_n = (2 * j + 1) / (4 * l + 2) * np.array(f_n)
                    self.channels.append(DiracChannel(k, f_n, basis))

        self.log('Basis functions: %s (%s)' %
                 (', '.join([str(nb) for nb in nb_l]),
                  ', '.join('spdf'[:lmax + 1])))

        self.vr_sg = self.rgd.zeros(self.nspins)
        self.vr_sg[:] = -self.Z

    def solve(self):
        """Diagonalize Schrödinger equation."""
        self.eeig = 0.0
        for channel in self.channels:
            if self.method == 'Gaussian basis-set':
                channel.solve(self.vr_sg[channel.s])
            else:
                channel.solve2(self.vr_sg[channel.s], self.scalar_relativistic,
                               self.Z)
            self.eeig += channel.get_eigenvalue_sum()

    def calculate_density(self):
        """Calculate elctron density and kinetic energy."""
        self.n_sg = self.rgd.zeros(self.nspins)
        for channel in self.channels:
            self.n_sg[channel.s] += channel.calculate_density()

    def calculate_electrostatic_potential(self):
        """Calculate electrostatic potential and energy."""
        n_g = self.n_sg.sum(0)
        self.vHr_g = self.rgd.poisson(n_g)
        self.eH = 0.5 * self.rgd.integrate(n_g * self.vHr_g, -1)
        self.eZ = -self.Z * self.rgd.integrate(n_g, -1)
        
    def calculate_xc_potential(self):
        self.vxc_sg = self.rgd.zeros(self.nspins)
        self.exc = self.xc.calculate_spherical(self.rgd, self.n_sg,
                                               self.vxc_sg)

    def step(self):
        self.solve()
        self.calculate_density()
        self.calculate_electrostatic_potential()
        self.calculate_xc_potential()
        self.vr_sg = self.vxc_sg * self.rgd.r_g
        self.vr_sg += self.vHr_g
        self.vr_sg -= self.Z
        self.ekin = (self.eeig -
                     self.rgd.integrate((self.vr_sg * self.n_sg).sum(0), -1))
        
    def run(self, mix=0.4, maxiter=117, dnmax=1e-9):
        if self.channels is None:
            self.initialize()

        if self.dirac:
            equation = 'Dirac'
        elif self.scalar_relativistic:
            equation = 'scalar-relativistic Schrödinger'
        else:
            equation = 'non-relativistic Schrödinger'
        self.log('\nSolving %s equation using %s:' % (equation, self.method))

        dn = self.Z

        vr_old_sg = None
        n_old_sg = None
        for iter in range(maxiter):
            self.log('.', end='')
            self.fd.flush()
            if iter > 0:
                self.vr_sg *= mix
                self.vr_sg += (1 - mix) * vr_old_sg
                dn = self.rgd.integrate(abs(self.n_sg - n_old_sg).sum(0))
                if dn <= dnmax:
                    self.log('\nConverged in', iter, 'steps')
                    break

            vr_old_sg = self.vr_sg
            n_old_sg = self.n_sg
            self.step()

        self.summary()
        
        if self.method != 'Gaussian basis-set':
            for channel in self.channels:
                assert channel.solve2ok

        if dn > dnmax:
            raise RuntimeError('Did not converge!')

    def refine(self):
        self.method = 'finite difference'
        self.run(dnmax=1e-6, mix=0.14, maxiter=200)
        
    def summary(self):
        self.write_states()
        self.write_energies()
            
    def write_states(self):
        self.log('\n state  occupation         eigenvalue          <r>')
        if self.dirac:
            self.log(' nl(j)               [Hartree]        [eV]    [Bohr]')
        else:
            self.log(' nl                  [Hartree]        [eV]    [Bohr]')
        self.log('-----------------------------------------------------')
        states = []
        for ch in self.channels:
            for n, f in enumerate(ch.f_n):
                states.append((ch.e_n[n], ch, n))
        states.sort()
        for e, ch, n in states:
            name = str(n + ch.l + 1) + ch.name
            if self.nspins == 2:
                name += '(%s)' % '+-'[ch.s]
            n_g = ch.calculate_density(n)
            rave = self.rgd.integrate(n_g, 1)
            self.log(' %-7s  %6.3f %13.6f  %13.5f %6.3f' %
                     (name, ch.f_n[n], e, e * units.Hartree, rave))

    def write_energies(self):
        self.log('\nEnergies:          [Hartree]           [eV]')
        self.log('--------------------------------------------')
        for text, e in [('kinetic      ', self.ekin),
                        ('coulomb (e-e)', self.eH),
                        ('coulomb (e-n)', self.eZ),
                        ('xc           ', self.exc),
                        ('total        ',
                         self.ekin + self.eH + self.eZ + self.exc)]:
            self.log(' %s %+13.6f  %+13.5f' % (text, e, units.Hartree * e))

        self.calculate_exx()
        self.log('\nExact exchange energy: %.6f Hartree, %.5f eV' %
                 (self.exx, self.exx * units.Hartree))

    def get_channel(self, l=None, s=0, k=None):
        if self.dirac:
            for channel in self.channels:
                if channel.k == k:
                    return channel
        else:
            for channel in self.channels:
                if channel.l == l and channel.s == s:
                    return channel
        raise ValueError

    def get_orbital(self, n, l=None, s=0, k=None):
        channel = self.get_channel(l, s, k)
        return channel.basis.expand(channel.C_nb[n])

    def plot_wave_functions(self, rc=4.0):
        import matplotlib.pyplot as plt
        for ch in self.channels:
            for n in range(len(ch.f_n)):
                fr_g = ch.basis.expand(ch.C_nb[n]) * self.rgd.r_g
                name = str(n + ch.l + 1) + ch.name
                lw = 2
                if self.nspins == 2:
                    name += '(%s)' % '+-'[ch.s]
                    if ch.s == 1:
                        lw = 1
                if self.dirac and ch.k > 0:
                    lw = 1
                ls = ['-', '--', '-.', ':'][ch.l]
                n_g = ch.calculate_density(n)
                rave = self.rgd.integrate(n_g, 1)
                gave = self.rgd.round(rave)
                fr_g *= cmp(fr_g[gave], 0)
                plt.plot(self.rgd.r_g, fr_g,
                         ls=ls, lw=lw, color=colors[n + ch.l], label=name)

        plt.legend(loc='best')
        plt.xlabel('r [Bohr]')
        plt.ylabel('$r\\phi(r)$')
        plt.axis(xmax=rc)
        plt.show()

    def logarithmic_derivative(self, l, energies, rcut):
        ch = Channel(l)
        gcut = self.rgd.round(rcut)
        u_g = self.rgd.empty()
        logderivs = []
        for e in energies:
            dudr = ch.integrate_outwards(u_g, self.rgd, self.vr_sg[0],
                                         gcut, e, self.scalar_relativistic)[0]
            logderivs.append(dudr / u_g[gcut])
        return logderivs
            
    def calculate_exx(self, s=None):
        if s is None:
            self.exx = sum(self.calculate_exx(s)
                           for s in range(self.nspins)) / self.nspins
            return self.exx

        states = []
        lmax = 0
        for ch in self.channels:
            l = ch.l
            for n, phi_g in enumerate(ch.phi_ng):
                f = ch.f_n[n]
                if f > 0 and ch.s == s:
                    states.append((l, f * self.nspins / 2.0 / (2 * l + 1),
                                   phi_g))
                    if l > lmax:
                        lmax = l

        G_LLL = make_gaunt(lmax)

        exx = 0.0
        j1 = 0
        for l1, f1, phi1_g in states:
            f = 1.0
            for l2, f2, phi2_g in states[j1:]:
                n_g = phi1_g * phi2_g
                for l in range((l1 + l2) % 2, l1 + l2 + 1, 2):
                    G = (G_LLL[l1**2:(l1 + 1)**2,
                               l2**2:(l2 + 1)**2,
                               l**2:(l + 1)**2]**2).sum()
                    vr_g = self.rgd.poisson(n_g, l)
                    e = f * self.rgd.integrate(vr_g * n_g, -1) / 4 / pi
                    exx -= e * G * f1 * f2
                f = 2.0
            j1 += 1

        return exx


def build_parser():
    from optparse import OptionParser

    parser = OptionParser(usage='gwap atom [options] element')
    parser.add_option('-f', '--xc-functional', type='string', default='LDA',
                      help='Exchange-Correlation functional ' +
                      '(default value LDA)',
                      metavar='<XC>')
    parser.add_option('-a', '--add', metavar='states',
                      help='Add electron(s). Use "1s0.5a" to add 0.5 1s ' +
                      'electrons to the alpha-spin channel (use "b" for ' +
                      'beta-spin).  The number of electrons defaults to ' +
                      'one. Examples: "1s", "2p2b", "4f0.1b,3d-0.1a".')
    parser.add_option('--spin-polarized', action='store_true')
    parser.add_option('-d', '--dirac', action='store_true')
    parser.add_option('-p', '--plot', action='store_true')
    parser.add_option('-e', '--exponents',
                      help='Exponents a: exp(-a*r^2).  Use "-e 0.1:20.0:30" ' +
                      'to get 30 exponents from 0.1 to 20.0.')
    parser.add_option('-l', '--logarithmic-derivatives',
                      metavar='spdfg,e1:e2:de,radius',
                      help='Plot logarithmic derivatives. ' +
                      'Example: -l spdf,-1:1:0.05,1.3. ' +
                      'Energy range and/or radius can be left out.')
    parser.add_option('-r', '--refine', action='store_true')
    parser.add_option('-s', '--scalar-relativistic', action='store_true')
    return parser


def parse_ld_str(s, energies=None, r=2.0):
    parts = s.split(',')
    lvalues = ['spdfg'.find(x) for x in parts.pop(0)]
    if parts:
        e1, e2, de = (float(x) for x in parts.pop(0).split(':'))
    else:
        e1, e2, de = energies
    if parts:
        r = float(parts.pop())
    energies = np.linspace(e1, e2, int((e2 - e1) / de) + 1)
    return lvalues, energies, r


def main(args=None):
    parser = build_parser()
    opt, args = parser.parse_args(args)

    if len(args) != 1:
        parser.error('Incorrect number of arguments')
    symbol = args[0]

    nlfs = []
    if opt.add:
        for x in opt.add.split(','):
            n = int(x[0])
            l = 'spdfg'.find(x[1])
            x = x[2:]
            if x and x[-1] in 'ab':
                s = int(x[-1] == 'b')
                opt.spin_polarized = True
                x = x[:-1]
            else:
                s = None
            if x:
                f = float(x)
            else:
                f = 1
            nlfs.append((n, l, f, s))

    aea = AllElectronAtom(symbol,
                          xc=opt.xc_functional,
                          spinpol=opt.spin_polarized,
                          dirac=opt.dirac)

    kwargs = {}
    if opt.exponents:
        parts = opt.exponents.split(':')
        kwargs['alpha1'] = float(parts[0])
        if len(parts) > 1:
            kwargs['alpha2'] = float(parts[1])
            if len(parts) > 2:
                kwargs['ngauss'] = int(parts[2])

    for n, l, f, s in nlfs:
        aea.add(n, l, f, s)

    aea.initialize(**kwargs)
    aea.run()

    if opt.refine:
        aea.refine()

    if opt.scalar_relativistic:
        aea.scalar_relativistic = True
        aea.refine()

    if opt.logarithmic_derivatives:
        lvalues, energies, r = parse_ld_str(opt.logarithmic_derivatives,
                                            (-1, 1, 0.05))
        import matplotlib.pyplot as plt
        for l in lvalues:
            ld = aea.logarithmic_derivative(l, energies, r)
            plt.plot(energies, ld, colors[l])
        plt.show()
        
    if opt.plot:
        aea.plot_wave_functions()


if __name__ == '__main__':
    main()
