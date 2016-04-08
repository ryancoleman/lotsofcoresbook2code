# -*- coding: utf-8 -*-
from __future__ import print_function
import sys
from math import pi, exp, sqrt, log

import numpy as np
from scipy.optimize import fsolve, root
from scipy import __version__ as scipy_version
from ase.units import Hartree
from ase.data import atomic_numbers

from gpaw.mpi import world
from gpaw.version import version
from gpaw.basis_data import Basis, BasisFunction, BasisPlotter
from gpaw.gaunt import make_gaunt
from gpaw.utilities import erf, pack2
from gpaw.utilities.lapack import general_diagonalize
from gpaw.atom.aeatom import AllElectronAtom, Channel, parse_ld_str, colors, \
    GaussianBasis


parameters = {
    # 1-2:
    'H1': ('1s,s,p', 0.9),
    'He2': ('1s,s,p', 1.5),
    # 3-10:
    'Li1': ('2s,s,2p', 2.1),
    'Li3': ('1s,2s,2p,p,d', 1.5),
    'Be2': ('2s,s,2p', 1.5),
    'Be4': ('1s,2s,2p,p,d', 1.4),
    'B3': ('2s,s,2p,p,d', 1.2),
    'C4': ('2s,s,2p,p,d', 1.2),
    'N5': ('2s,s,2p,p,d', [1.2, 1.3], {'r0': 1.1}),
    'O6': ('2s,s,2p,p,d,F', 1.2),
    'F7': ('2s,s,2p,p,d', [1.2, 1.4]),
    'Ne8': ('2s,s,2p,p,d', 1.8),
    # 11-18:
    'Na1': ('3s,s,3p', 2.6),
    'Na9': ('2s,3s,2p,3p,d,F', 2.3),
    'Mg2': ('3s,s,3p,D', 2.6),
    'Mg10': ('2s,3s,2p,3p,d,F', [2.0, 1.8]),
    'Al3': ('3s,s,3p,p,d,F', 2.1),
    'Si4': ('3s,s,3p,p,d,F', 1.9),
    'P5': ('3s,s,3p,p,d,F', 1.7),
    'S6': ('3s,s,3p,p,d,F', 1.6),
    'Cl7': ('3s,s,3p,p,d,F', 1.5),
    'Ar8': ('3s,s,3p,p,d,F', 1.5),
    # 19-36:
    'K1': ('4s,s,4p,D', 3.5),
    'K9': ('3s,4s,3p,4p,d,d,F', 2.1),
    'Ca2': ('4s,s,4p', 3.1),
    'Ca10': ('3s,4s,3p,4p,3d,d,F', 2.1),
    'Sc3': ('4s,s,4p,p,3d,d', 2.7),
    'Sc11': ('3s,4s,3p,4p,3d,d,F', 2.3),
    'Ti4': ('4s,s,4p,p,3d,d', 2.7),
    'Ti12': ('3s,4s,3p,4p,3d,d,F', [2.2, 2.2, 2.3]),
    'V5': ('4s,s,4p,p,3d,d', 2.6),
    'V13': ('3s,4s,3p,4p,3d,d,F', [2.1, 2.1, 2.3]),
    'Cr6': ('4s,s,4p,p,3d,d', 2.5),
    'Cr14': ('3s,4s,3p,4p,3d,d,F', [2.1, 2.1, 2.3]),
    'Mn7': ('4s,s,4p,p,3d,d', 2.4),
    'Mn15': ('3s,4s,3p,4p,3d,d,F', [2.0, 2.0, 2.2]),
    'Fe8': ('4s,s,4p,p,3d,d', 2.2),
    'Fe16': ('3s,4s,3p,4p,3d,d,F', 2.1),
    'Co9': ('4s,s,4p,p,3d,d', 2.2),
    'Co17': ('3s,4s,3p,4p,3d,d,F', 2.1),
    'Ni10': ('4s,s,4p,p,3d,d', 2.1),
    'Ni18': ('3s,4s,3p,4p,3d,d,F', 2.0),
    'Cu11': ('4s,s,4p,p,3d,d', 2.1),
    'Cu19': ('3s,4s,3p,4p,3d,d,F', 1.9),
    'Zn12': ('4s,s,4p,p,3d', 2.1),
    'Zn20': ('3s,4s,3p,4p,3d,d,F', 1.9),
    'Ga3': ('4s,s,4p,p,d,F', 2.2),
    'Ga13': ('4s,s,4p,p,3d,d,F', 2.2),
    'Ge4': ('4s,s,4p,p,d,F', 2.1),
    'Ge14': ('4s,s,4p,p,3d,d,F', 2.1),
    'As5': ('4s,s,4p,p,d,F', 2.0),
    'Se6': ('4s,s,4p,p,d,F', 2.1),
    'Br7': ('4s,s,4p,p,d,F', 2.1),
    'Kr8': ('4s,s,4p,p,d,F', 2.1),
    # 37-54:
    'Rb1': ('5s,s,5p', 3.6),
    'Rb9': ('4s,5s,4p,5p,d,d,F', 2.5),
    'Sr2': ('5s,s,5p', 3.3),
    'Sr10': ('4s,5s,4p,5p,4d,d,F', 2.5),
    'Y3': ('5s,s,5p,p,4d,d', 3.1),
    'Y11': ('4s,5s,4p,5p,4d,d,F', 2.5),
    'Zr4': ('5s,s,5p,p,4d,d', 3.0),
    'Zr12': ('4s,5s,4p,5p,4d,d,F', 2.5),
    'Nb5': ('5s,s,5p,p,4d,d', 2.9),
    'Nb13': ('4s,5s,4p,5p,4d,d,F', [2.4, 2.4, 2.5]),
    'Mo6': ('5s,s,5p,p,4d,d', 2.8),
    'Mo14': ('4s,5s,4p,5p,4d,d,F', 2.3),
    'Tc7': ('5s,s,5p,p,4d,d', 2.7),
    'Tc15': ('4s,5s,4p,5p,4d,d,F', 2.3),
    'Ru8': ('5s,s,5p,p,4d,d', 2.6),
    'Ru16': ('4s,5s,4p,5p,4d,d,F', 2.3),
    'Rh9': ('5s,s,5p,p,4d,d', 2.5),
    'Rh17': ('4s,5s,4p,5p,4d,d,F', 2.3),
    'Pd10': ('5s,s,5p,p,4d,d', 2.4),
    'Pd18': ('4s,5s,4p,5p,4d,d,F', 2.3),
    'Ag11': ('5s,s,5p,p,4d,d', 2.4),
    'Ag19': ('4s,5s,4p,5p,4d,d,F', 2.3),
    'Cd12': ('5s,s,5p,p,4d,d', 2.4),
    'Cd20': ('4s,5s,4p,5p,4d,d,F', 2.3),
    'In13': ('5s,s,5p,p,4d,d,F', 2.6),
    'Sn14': ('5s,s,5p,p,4d,d,F', 2.5),
    'Sb15': ('5s,s,5p,p,4d,d,F', 2.5),
    'Te6': ('5s,6s,5p,p,d,d,F', 2.5),
    'I7': ('5s,s,5p,p,d,F', 2.4),
    'Xe8': ('5s,s,5p,p,d,F', 2.3),
    # 55-56:
    'Cs1': ('6s,s,6p,5d', [4.3, 4.6, 4.0]),
    'Cs9': ('5s,6s,5p,6p,5d', [1.9, 2.2]),
    'Ba2': ('6s,s,6p,5d', 3.9),
    'Ba10': ('5s,6s,5p,6p,5d', [1.8, 2.2]),
    # 57-71:
    'La11': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.5),
    'Ce12': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.4),
    'Pr13': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.3),
    'Nd14': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.3),
    'Pm15': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.3),
    'Sm16': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.2),
    'Eu17': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.2),
    'Gd18': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.2),
    'Tb19': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.2),
    'Dy20': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.1),
    'Ho21': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.2),
    'Er22': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.2),
    'Tm23': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.2),
    'Yb24': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.2),
    'Lu25': ('5s,6s,5p,6p,5d,d,4f,f,G', 2.2),
    # 72-86:
    'Hf4': ('6s,s,6p,p,5d,d', 2.9),
    'Hf12': ('5s,6s,5p,6p,5d,d,F', 2.4),
    'Ta5': ('6s,s,6p,p,5d,d', 2.8),
    'Ta13': ('5s,6s,5p,6p,5d,d,F', 2.4),
    'W6': ('6s,s,6p,p,5d,d', 2.7),
    'W14': ('5s,6s,5p,6p,5d,d,F', 2.4),
    'Re7': ('6s,s,6p,p,5d,d', 2.6),
    'Re15': ('5s,6s,5p,6p,5d,d,F', 2.4),
    'Os8': ('6s,s,6p,p,5d,d', 2.6),
    'Os16': ('5s,6s,5p,6p,5d,d,F', 2.4),
    'Ir9': ('6s,s,6p,p,5d,d', 2.6),
    'Ir17': ('5s,6s,5p,6p,5d,d,F', 2.4),
    'Pt10': ('6s,s,6p,p,5d,d', 2.5),
    'Pt18': ('5s,6s,5p,6p,5d,d,F', 2.3),
    'Au11': ('6s,s,6p,p,5d,d', 2.5),
    'Au19': ('5s,6s,5p,6p,5d,d,F', 2.3),
    'Hg12': ('6s,s,6p,p,5d,d', 2.5),
    'Hg20': ('5s,6s,5p,6p,5d,d,F', 2.3),
    'Tl13': ('6s,s,6p,p,5d,d,F', 2.8),
    'Pb14': ('6s,s,6p,p,5d,d,F', 2.6),
    'Bi5': ('6s,s,6p,p,d,F', 2.8),
    'Bi15': ('6s,s,6p,p,5d,d,F', 2.6),
    'Po6': ('6s,s,6p,p,d,F', 2.7),
    'At7': ('6s,s,6p,p,d,F', 2.6),
    'Rn8': ('6s,s,6p,p,d,F', 2.6)}


default = [0,
           1, 2,
           1, 2, 3, 4, 5, 6, 7, 8,
           1, 2, 3, 4, 5, 6, 7, 8,
           1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8,
           1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 6, 7, 8,
           1, 2, 11,
           12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
           4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 6, 7, 8]


semicore = [0,
            1, 2,
            3, 4, 3, 4, 5, 6, 7, 8,
            9, 10, 3, 4, 5, 6, 7, 8,
            9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 13, 14, 5, 6, 7, 8,
            9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 13, 14, 15, 6, 7, 8,
            9, 10, 11,
            12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
            12, 13, 14, 15, 16, 17, 18, 19, 20, 13, 14, 15, 6, 7, 8]


def get_number_of_electrons(symbol, name):
    Z = atomic_numbers[symbol]
    if name == 'default':
        return default[Z]
    return semicore[Z]
    
        
class PAWWaves:
    def __init__(self, rgd, l, rcut):
        self.rgd = rgd
        self.l = l
        self.rcut = rcut

        self.n_n = []
        self.e_n = []
        self.f_n = []
        self.phi_ng = []
        self.phit_ng = None
        self.pt_ng = None

    def __len__(self):
        return len(self.n_n)

    def add(self, phi_g, n, e, f):
        self.phi_ng.append(phi_g)
        self.n_n.append(n)
        self.e_n.append(e)
        self.f_n.append(f)

    def pseudize(self, type, nderiv):
        rgd = self.rgd

        if type == 'poly':
            ps = rgd.pseudize
        elif type == 'bessel':
            ps = rgd.jpseudize

        phi_ng = self.phi_ng = np.array(self.phi_ng)
        N = len(phi_ng)
        phit_ng = self.phit_ng = rgd.empty(N)
        gcut = rgd.ceil(self.rcut)

        self.nt_g = 0
        self.c_n = []
        for n in range(N):
            phit_ng[n], c = ps(phi_ng[n], gcut, self.l, nderiv)
            self.c_n.append(c)
            self.nt_g += self.f_n[n] / 4 / pi * phit_ng[n]**2

        self.dS_nn = np.empty((N, N))
        for n1 in range(N):
            for n2 in range(N):
                self.dS_nn[n1, n2] = rgd.integrate(
                    phi_ng[n1] * phi_ng[n2] -
                    phit_ng[n1] * phit_ng[n2]) / (4 * pi)
        self.Q = np.dot(self.f_n, self.dS_nn.diagonal())

    def construct_projectors(self, vtr_g, rcmax):
        N = len(self)
        if N == 0:
            self.pt_ng = []
            return

        rgd = self.rgd
        phit_ng = self.phit_ng
        gcmax = rgd.ceil(rcmax)
        r_g = rgd.r_g
        l = self.l

        dgdr_g = 1 / rgd.dr_g
        d2gdr2_g = rgd.d2gdr2()

        q_ng = rgd.zeros(N)
        for n in range(N):
            a_g, dadg_g, d2adg2_g = rgd.zeros(3)
            a_g[1:] = self.phit_ng[n, 1:] / r_g[1:]**l
            a_g[0] = self.c_n[n]
            dadg_g[1:-1] = 0.5 * (a_g[2:] - a_g[:-2])
            d2adg2_g[1:-1] = a_g[2:] - 2 * a_g[1:-1] + a_g[:-2]
            q_g = (vtr_g - self.e_n[n] * r_g) * self.phit_ng[n]
            q_g -= 0.5 * r_g**l * (
                (2 * (l + 1) * dgdr_g + r_g * d2gdr2_g) * dadg_g +
                r_g * d2adg2_g * dgdr_g**2)
            q_g[gcmax:] = 0
            q_g[1:] /= r_g[1:]
            if l == 0:
                q_g[0] = q_g[1]
            q_ng[n] = q_g

        A_nn = rgd.integrate(phit_ng[:, None] * q_ng) / (4 * pi)
        self.dH_nn = self.e_n * self.dS_nn - A_nn

        L_nn = np.eye(N)
        U_nn = A_nn.copy()

        if N - self.n_n.count(-1) == 1:
            assert self.n_n[0] != -1
            # We have a single bound-state projector.
            for n1 in range(N):
                for n2 in range(n1 + 1, N):
                    L_nn[n2, n1] = U_nn[n2, n1] / U_nn[n1, n1]
                    U_nn[n2] -= U_nn[n1] * L_nn[n2, n1]

            iL_nn = np.linalg.inv(L_nn)
            phit_ng[:] = np.dot(iL_nn, phit_ng)
            self.phi_ng[:] = np.dot(iL_nn, self.phi_ng)

            self.dS_nn = np.dot(np.dot(iL_nn, self.dS_nn), iL_nn.T)
            self.dH_nn = np.dot(np.dot(iL_nn, self.dH_nn), iL_nn.T)

        self.pt_ng = np.dot(np.linalg.inv(U_nn.T), q_ng)

    def calculate_kinetic_energy_correction(self, vr_g, vtr_g):
        if len(self) == 0:
            return
        self.dekin_nn = (self.rgd.integrate(self.phit_ng[:, None] *
                                            self.phit_ng *
                                            vtr_g, -1) / (4 * pi) -
                         self.rgd.integrate(self.phi_ng[:, None] *
                                            self.phi_ng *
                                            vr_g, -1) / (4 * pi) +
                         self.dH_nn)


class PAWSetupGenerator:
    def __init__(self, aea, projectors,
                 scalar_relativistic=False,
                 core_hole=None,
                 fd=None):
        """fd: stream
            Text output."""

        self.aea = aea

        self.fd = fd or sys.stdout

        if core_hole:
            state, occ = core_hole.split(',')
            n = int(state[0])
            l = 'spdf'.find(state[1])
            occ = float(occ)
            aea.add(n, l, -occ)
            self.core_hole = (n, l, occ)
        else:
            self.core_hole = None

        if projectors[-1].isupper():
            self.l0 = 'SPDFG'.find(projectors[-1])
            projectors = projectors[:-2]
        else:
            self.l0 = None

        self.lmax = -1
        self.states = {}
        for s in projectors.split(','):
            l = 'spdf'.find(s[-1])
            if len(s) == 1:
                n = None
            elif '.' in s:
                n = float(s[:-1])
            else:
                n = int(s[:-1])
            if l in self.states:
                self.states[l].append(n)
            else:
                self.states[l] = [n]
            if l > self.lmax:
                self.lmax = l

        # Add empty bound states:
        for l, nn in self.states.items():
            for n in nn:
                if (isinstance(n, int) and
                    (l not in aea.f_lsn or n - l > len(aea.f_lsn[l][0]))):
                    aea.add(n, l, 0)

        aea.initialize()
        aea.run()
        aea.scalar_relativistic = scalar_relativistic
        aea.refine()

        self.rgd = aea.rgd

        self.vtr_g = None

        self.basis = None

        self.log('\nGenerating PAW', aea.xc.name, 'setup for', aea.symbol)

    def construct_shape_function(self, alpha=None, rc=None, eps=1e-10):
        """Build shape-function for compensation charge."""

        self.alpha = alpha

        if self.alpha is None:
            if isinstance(rc, list):
                rc = min(rc)
            rc = 1.5 * rc

            def spillage(alpha):
                """Fraction of gaussian charge outside rc."""
                x = alpha * rc**2
                return 1 - erf(sqrt(x)) + 2 * sqrt(x / pi) * exp(-x)

            def f(alpha):
                return log(spillage(alpha)) - log(eps)

            if scipy_version < '0.8':
                self.alpha = fsolve(f, 7.0)
            else:
                self.alpha = fsolve(f, 7.0)[0]

            self.alpha = round(self.alpha, 1)

        self.log('Shape function: exp(-alpha*r^2), alpha=%.1f Bohr^-2' %
                 self.alpha)

        self.ghat_g = (np.exp(-self.alpha * self.rgd.r_g**2) *
                       (self.alpha / pi)**1.5)

    def log(self, *args, **kwargs):
        print(file=self.fd, *args, **kwargs)

    def calculate_core_density(self):
        self.nc_g = self.rgd.zeros()
        self.tauc_g = self.rgd.zeros()
        self.ncore = 0
        self.nvalence = 0
        self.ekincore = 0.0
        for l, ch in enumerate(self.aea.channels):
            for n, f in enumerate(ch.f_n):
                if (l <= self.lmax and
                    any(n + l + 1 == nn
                        for nn in self.states[l]
                        if isinstance(nn, int))):
                    self.nvalence += f
                else:
                    self.nc_g += f * ch.calculate_density(n)
                    self.tauc_g += f * ch.calculate_kinetic_energy_density(n)
                    self.ncore += f
                    self.ekincore += f * ch.e_n[n]

        self.ekincore -= self.rgd.integrate(self.nc_g * self.aea.vr_sg[0], -1)

        self.log('Core electrons:', self.ncore)
        self.log('Valence electrons:', self.nvalence)

    def add_waves(self, rc):
        if isinstance(rc, float):
            radii = [rc]
        else:
            radii = rc

        self.rcmax = max(radii)

        if self.lmax >= 0:
            radii += [radii[-1]] * (self.lmax + 1 - len(radii))

        self.waves_l = []
        for l in range(self.lmax + 1):
            rcut = radii[l]
            waves = PAWWaves(self.rgd, l, rcut)
            e = -1.0
            for n in self.states[l]:
                if isinstance(n, int):
                    # Bound state:
                    ch = self.aea.channels[l]
                    e = ch.e_n[n - l - 1]
                    f = ch.f_n[n - l - 1]
                    phi_g = ch.phi_ng[n - l - 1]
                else:
                    if n is None:
                        e += 1.0
                    else:
                        e = n
                    n = -1
                    f = 0.0
                    phi_g = self.rgd.zeros()
                    gc = self.rgd.round(1.5 * rcut)
                    ch = Channel(l)
                    a = ch.integrate_outwards(phi_g, self.rgd,
                                              self.aea.vr_sg[0], gc, e,
                                              self.aea.scalar_relativistic,
                                              self.aea.Z)[1]
                    phi_g[1:gc + 1] /= self.rgd.r_g[1:gc + 1]
                    phi_g[0] = a
                    phi_g /= (self.rgd.integrate(phi_g**2) / (4 * pi))**0.5

                waves.add(phi_g, n, e, f)
            self.waves_l.append(waves)

    def pseudize(self, type='poly', nderiv=6, rcore=None):
        self.Q = -self.aea.Z + self.ncore

        self.nt_g = self.rgd.zeros()
        for waves in self.waves_l:
            waves.pseudize(type, nderiv)
            self.nt_g += waves.nt_g
            self.Q += waves.Q

        if rcore is None:
            rcore = self.rcmax * 0.8
        else:
            assert rcore <= self.rcmax

        if self.ncore == 0:
            self.nct_g = self.rgd.zeros()
            self.tauct_g = self.rgd.zeros()
        else:
            # Make sure pseudo density is monotonically decreasing:
            while 1:
                gcore = self.rgd.round(rcore)
                self.nct_g = self.rgd.pseudize(self.nc_g, gcore)[0]
                nt_g = self.nt_g + self.nct_g
                dntdr_g = self.rgd.derivative(nt_g)[:gcore]
                if dntdr_g.max() < 0.0:
                    break
                rcore -= 0.01

            rcore *= 1.2
            gcore = self.rgd.round(rcore)
            self.nct_g = self.rgd.pseudize(self.nc_g, gcore)[0]
            nt_g = self.nt_g + self.nct_g

            self.log('Constructing smooth pseudo core density for r < %.3f' %
                     rcore)
            self.nt_g = nt_g

            self.tauct_g = self.rgd.pseudize(self.tauc_g, gcore)[0]

        self.npseudocore = self.rgd.integrate(self.nct_g)
        self.log('Pseudo core electrons: %.6f' % self.npseudocore)
        self.Q -= self.npseudocore

        self.rhot_g = self.nt_g + self.Q * self.ghat_g
        self.vHtr_g = self.rgd.poisson(self.rhot_g)

        self.vxct_g = self.rgd.zeros()
        nt_sg = self.nt_g.reshape((1, -1))
        self.exct = self.aea.xc.calculate_spherical(
            self.rgd, nt_sg, self.vxct_g.reshape((1, -1)))

        self.v0r_g = self.vtr_g - self.vHtr_g - self.vxct_g * self.rgd.r_g
        self.v0r_g[self.rgd.round(self.rcmax):] = 0.0

        self.log('\nProjectors:')
        self.log(' state  occ         energy             norm        rcut')
        self.log(' nl            [Hartree]  [eV]      [electrons]   [Bohr]')
        self.log('----------------------------------------------------------')
        for l, waves in enumerate(self.waves_l):
            for n, e, f, ds in zip(waves.n_n, waves.e_n, waves.f_n,
                                   waves.dS_nn.diagonal()):
                if n == -1:
                    self.log('  %s         %10.6f %10.5f   %19.2f' %
                             ('spdf'[l], e, e * Hartree, waves.rcut))
                else:
                    self.log(
                        ' %d%s   %5.2f %10.6f %10.5f      %5.3f  %9.2f' %
                        (n, 'spdf'[l], f, e, e * Hartree, 1 - ds,
                         waves.rcut))
        self.log()

    def find_local_potential(self, r0, P):
        self.r0 = r0
        self.nderiv0 = P
        if self.l0 is None:
            self.find_polynomial_potential(r0, P)
        else:
            self.match_local_potential(r0, P)

    def find_polynomial_potential(self, r0, P):
        self.log('Constructing smooth local potential for r < %.3f' % r0)
        g0 = self.rgd.ceil(r0)

        r_g = self.rgd.r_g
        if 1:
            self.vtr_g = self.rgd.pseudize(self.aea.vr_sg[0], g0, 1, P)[0]
        else:
            # Use bessel function:
            r0 = r_g[g0]
            v_g = self.aea.vr_sg[0, g0 - 1:g0 + 2] / r_g[g0 - 1:g0 + 2]
            v0 = v_g[1]
            v1 = (v_g[2] - v_g[0]) / 2 / self.rgd.dr_g[g0]

            def f(x):
                return x / np.tan(x) - 1 - v1 * r0 / v0
            
            q = root(f, 2).x / r0
            A = v0 * r0 / np.sin(q * r0)
            self.vtr_g = self.aea.vr_sg[0].copy()
            self.vtr_g[:g0] = A * np.sin(q * r_g[:g0])

    def match_local_potential(self, r0, P):
        l0 = self.l0
        self.log('Local potential matching %s-scattering at e=0.0 eV' %
                 'spdfg'[l0] + ' and r=%.2f Bohr' % r0)

        g0 = self.rgd.ceil(r0)
        gc = g0 + 20

        e0 = 0.0

        ch = Channel(l0)
        phi_g = self.rgd.zeros()
        a = ch.integrate_outwards(phi_g, self.rgd, self.aea.vr_sg[0], gc, e0,
                                  self.aea.scalar_relativistic, self.aea.Z)[1]
        phi_g[1:gc] /= self.rgd.r_g[1:gc]
        phi_g[0] = a

        phit_g, c = self.rgd.pseudize(phi_g, g0, l=l0, points=P)

        dgdr_g = 1 / self.rgd.dr_g
        d2gdr2_g = self.rgd.d2gdr2()
        a_g = phit_g.copy()
        a_g[1:] /= self.rgd.r_g[1:]**l0
        a_g[0] = c
        dadg_g = self.rgd.zeros()
        d2adg2_g = self.rgd.zeros()
        dadg_g[1:-1] = 0.5 * (a_g[2:] - a_g[:-2])
        d2adg2_g[1:-1] = a_g[2:] - 2 * a_g[1:-1] + a_g[:-2]
        q_g = (((l0 + 1) * dgdr_g + 0.5 * self.rgd.r_g * d2gdr2_g) * dadg_g +
               0.5 * self.rgd.r_g * d2adg2_g * dgdr_g**2)
        q_g[:g0] /= a_g[:g0]
        q_g += e0 * self.rgd.r_g
        q_g[0] = 0.0

        self.vtr_g = self.aea.vr_sg[0].copy()
        self.vtr_g[0] = 0.0
        self.vtr_g[1:g0] = q_g[1:g0]

    def construct_projectors(self):
        for waves in self.waves_l:
            waves.construct_projectors(self.vtr_g, self.rcmax)
            waves.calculate_kinetic_energy_correction(self.aea.vr_sg[0],
                                                      self.vtr_g)

    def check_all(self):
        self.log(('Checking eigenvalues of %s pseudo atom using ' +
                  'a Gaussian basis set:') % self.aea.symbol)
        self.log('                 AE [eV]        PS [eV]      error [eV]')

        ok = True

        for l in range(4):
            try:
                e_b, n0 = self.check(l)
            except RuntimeError:
                self.log('Singular overlap matrix!')
                ok = False
                continue

            nbound = (e_b < -0.002).sum()

            if l < len(self.aea.channels):
                e0_b = self.aea.channels[l].e_n
                nbound0 = (e0_b < -0.002).sum()
                extra = 6
                for n in range(1 + l, nbound0 + 1 + l + extra):
                    if n - 1 - l < len(self.aea.channels[l].f_n):
                        f = self.aea.channels[l].f_n[n - 1 - l]
                        self.log('%2d%s  %2d' % (n, 'spdf'[l], f), end='')
                    else:
                        self.log('       ', end='')
                    self.log('  %15.3f' % (e0_b[n - 1 - l] * Hartree), end='')
                    if n - 1 - l - n0 >= 0:
                        self.log('%15.3f' * 2 %
                                 (e_b[n - 1 - l - n0] * Hartree,
                                  (e_b[n - 1 - l - n0] - e0_b[n - 1 - l]) *
                                  Hartree))
                    else:
                        self.log()

                if nbound != nbound0 - n0:
                    self.log('Wrong number of %s-states!' % 'spdf'[l])
                    ok = False
                elif (nbound > 0 and
                      abs(e_b[:nbound] - e0_b[n0:nbound0]).max() > 1e-3):
                    self.log('Error in bound %s-states!' % 'spdf'[l])
                    ok = False
                elif (abs(e_b[nbound:nbound + extra] -
                          e0_b[nbound0:nbound0 + extra]).max() > 2e-2):
                    self.log('Error in %s-states!' % 'spdf'[l])
                    ok = False
            elif nbound > 0:
                self.log('Wrong number of %s-states!' % 'spdf'[l])
                ok = False

        return ok

    def check(self, l):
        basis = self.aea.channels[0].basis
        eps = basis.eps
        alpha_B = basis.alpha_B

        basis = GaussianBasis(l, alpha_B, self.rgd, eps)
        H_bb = basis.calculate_potential_matrix(self.vtr_g)
        H_bb += basis.T_bb
        S_bb = np.eye(len(basis))

        n0 = 0
        if l < len(self.waves_l):
            waves = self.waves_l[l]
            if len(waves) > 0:
                P_bn = self.rgd.integrate(basis.basis_bg[:, None] *
                                          waves.pt_ng) / (4 * pi)
                H_bb += np.dot(np.dot(P_bn, waves.dH_nn), P_bn.T)
                S_bb += np.dot(np.dot(P_bn, waves.dS_nn), P_bn.T)
                n0 = waves.n_n[0] - l - 1
                if n0 < 0 and l < len(self.aea.channels):
                    n0 = (self.aea.channels[l].f_n > 0).sum()
        elif l < len(self.aea.channels):
            n0 = (self.aea.channels[l].f_n > 0).sum()

        e_b = np.empty(len(basis))
        general_diagonalize(H_bb, e_b, S_bb)
        return e_b, n0

    def test_convergence(self):
        rgd = self.rgd
        r_g = rgd.r_g
        G_k, nt_k = self.rgd.fft(self.nt_g * r_g)
        rhot_k = self.rgd.fft(self.rhot_g * r_g)[1]
        ghat_k = self.rgd.fft(self.ghat_g * r_g)[1]
        vt_k = self.rgd.fft(self.vtr_g)[1]
        phi_k = self.rgd.fft(self.waves_l[0].phit_ng[0] * r_g)[1]
        eee_k = 0.5 * nt_k**2 * (4 * pi)**2 / (2 * pi)**3
        ecc_k = 0.5 * rhot_k**2 * (4 * pi)**2 / (2 * pi)**3
        egg_k = 0.5 * ghat_k**2 * (4 * pi)**2 / (2 * pi)**3
        ekin_k = 0.5 * phi_k**2 * G_k**4 / (2 * pi)**3
        evt_k = nt_k * vt_k * G_k**2 * 4 * pi / (2 * pi)**3

        eee = 0.5 * rgd.integrate(self.nt_g * rgd.poisson(self.nt_g), -1)
        ecc = 0.5 * rgd.integrate(self.rhot_g * self.vHtr_g, -1)
        egg = 0.5 * rgd.integrate(self.ghat_g * rgd.poisson(self.ghat_g), -1)
        ekin = self.aea.ekin - self.ekincore - self.waves_l[0].dekin_nn[0, 0]
        evt = rgd.integrate(self.nt_g * self.vtr_g, -1)

        import pylab as p

        errors = 10.0**np.arange(-4, 0) / Hartree
        self.log('\nConvergence of energy:')
        self.log('plane-wave cutoff (wave-length) [ev (Bohr)]\n  ', end='')
        for de in errors:
            self.log('%14.4f' % (de * Hartree), end='')
        for label, e_k, e0 in [
            ('e-e', eee_k, eee),
            ('c-c', ecc_k, ecc),
            ('g-g', egg_k, egg),
            ('kin', ekin_k, ekin),
            ('vt', evt_k, evt)]:
            self.log('\n%3s: ' % label, end='')
            e_k = (np.add.accumulate(e_k) - 0.5 * e_k[0] - 0.5 * e_k) * G_k[1]
            k = len(e_k) - 1
            for de in errors:
                while abs(e_k[k] - e_k[-1]) < de:
                    k -= 1
                G = k * G_k[1]
                ecut = 0.5 * G**2
                h = pi / G
                self.log(' %6.1f (%4.2f)' % (ecut * Hartree, h), end='')
            p.semilogy(G_k, abs(e_k - e_k[-1]) * Hartree, label=label)
        self.log()
        p.axis(xmax=20)
        p.xlabel('G')
        p.ylabel('[eV]')
        p.legend()
        p.show()

    def plot(self):
        import matplotlib.pyplot as plt
        r_g = self.rgd.r_g

        plt.figure()
        plt.plot(r_g, self.vxct_g, label='xc')
        plt.plot(r_g[1:], self.v0r_g[1:] / r_g[1:], label='0')
        plt.plot(r_g[1:], self.vHtr_g[1:] / r_g[1:], label='H')
        plt.plot(r_g[1:], self.vtr_g[1:] / r_g[1:], label='ps')
        plt.plot(r_g[1:], self.aea.vr_sg[0, 1:] / r_g[1:], label='ae')
        plt.axis(xmax=2 * self.rcmax,
                 ymin=self.vtr_g[1] / r_g[1],
                 ymax=max(0, (self.v0r_g[1:] / r_g[1:]).max()))
        plt.xlabel('radius [Bohr]')
        plt.ylabel('potential [Ha]')
        plt.legend()

        plt.figure()
        i = 0
        for l, waves in enumerate(self.waves_l):
            for n, e, phi_g, phit_g in zip(waves.n_n, waves.e_n,
                                           waves.phi_ng, waves.phit_ng):
                if n == -1:
                    gc = self.rgd.ceil(waves.rcut)
                    name = '*%s (%.2f Ha)' % ('spdf'[l], e)
                else:
                    gc = len(self.rgd)
                    name = '%d%s (%.2f Ha)' % (n, 'spdf'[l], e)
                plt.plot(r_g[:gc], (phi_g * r_g)[:gc], color=colors[i],
                         label=name)
                plt.plot(r_g[:gc], (phit_g * r_g)[:gc], '--', color=colors[i])
                i += 1
        plt.axis(xmax=3 * self.rcmax)
        plt.xlabel('radius [Bohr]')
        plt.ylabel(r'$r\phi_{n\ell}(r)$')
        plt.legend()

        plt.figure()
        i = 0
        for l, waves in enumerate(self.waves_l):
            for n, e, pt_g in zip(waves.n_n, waves.e_n, waves.pt_ng):
                if n == -1:
                    name = '*%s (%.2f Ha)' % ('spdf'[l], e)
                else:
                    name = '%d%s (%.2f Ha)' % (n, 'spdf'[l], e)
                plt.plot(r_g, pt_g * r_g, color=colors[i], label=name)
                i += 1
        plt.axis(xmax=self.rcmax)
        plt.legend()

    def create_basis_set(self, tailnorm=0.0005, scale=200.0, splitnorm=0.16):
        rgd = self.rgd
        self.basis = Basis(self.aea.symbol, readxml=False, rgd=rgd)

        # We print text to sdtout and put it in the basis-set file
        txt = 'Basis functions:\n'

        # Bound states:
        for l, waves in enumerate(self.waves_l):
            for i, n in enumerate(waves.n_n):
                if n > 0:
                    tn = tailnorm
                    if waves.f_n[i] == 0:
                        tn = min(0.05, tn * 20)  # no need for long tail
                    phit_g, ronset, rc, de = self.create_basis_function(
                        l, i, tn, scale)
                    bf = BasisFunction(n, l, rc, phit_g, 'bound state')
                    self.basis.append(bf)

                    txt += '%d%s bound state:\n' % (n, 'spdf'[l])
                    txt += ('  cutoff: %.3f to %.3f Bohr (tail-norm=%f)\n' %
                            (ronset, rc, tn))
                    txt += '  eigenvalue shift: %.3f eV\n' % (de * Hartree)

        # Split valence:
        for l, waves in enumerate(self.waves_l):
            # Find the largest n that is occupied:
            n0 = None
            for f, n in zip(waves.f_n, waves.n_n):
                if n > 0 and f > 0:
                    n0 = n
            if n0 is None:
                continue

            bf = [bf for bf in self.basis.bf_j if bf.l == l and bf.n == n0][0]

            # Radius and l-value used for polarization function below:
            rcpol = bf.rc
            lpol = l + 1

            phit_g = bf.phit_g

            # Find cutoff radius:
            n_g = np.add.accumulate(phit_g**2 * rgd.r_g**2 * rgd.dr_g)
            norm = n_g[-1]
            gc = (norm - n_g > splitnorm * norm).sum()
            rc = rgd.r_g[gc]

            phit2_g = rgd.pseudize(phit_g, gc, l, 2)[0]  # "split valence"
            bf = BasisFunction(n, l, rc, phit_g - phit2_g, 'split valence')
            self.basis.append(bf)

            txt += '%d%s split valence:\n' % (n0, 'spdf'[l])
            txt += '  cutoff: %.3f Bohr (tail-norm=%f)\n' % (rc, splitnorm)

        # Polarization:
        gcpol = rgd.round(rcpol)
        alpha = 1 / (0.25 * rcpol)**2

        # Gaussian that is continuous and has a continuous derivative at rcpol:
        phit_g = np.exp(-alpha * rgd.r_g**2) * rgd.r_g**lpol
        phit_g -= rgd.pseudize(phit_g, gcpol, lpol, 2)[0]
        phit_g[gcpol:] = 0.0

        bf = BasisFunction(None, lpol, rcpol, phit_g, 'polarization')
        self.basis.append(bf)
        txt += 'l=%d polarization functions:\n' % lpol
        txt += '  cutoff: %.3f Bohr (r^%d exp(-%.3f*r^2))\n' % (rcpol, lpol,
                                                                alpha)

        self.log(txt)

        # Write basis-set file:
        self.basis.generatordata = txt
        self.basis.generatorattrs.update(dict(tailnorm=tailnorm,
                                              scale=scale,
                                              splitnorm=splitnorm))
        self.basis.name = '%de.dzp' % self.nvalence

        return self.basis

    def create_basis_function(self, l, n, tailnorm, scale):
        rgd = self.rgd
        waves = self.waves_l[l]

        # Find cutoff radii:
        n_g = np.add.accumulate(waves.phit_ng[n]**2 * rgd.r_g**2 * rgd.dr_g)
        norm = n_g[-1]
        g2 = (norm - n_g > tailnorm * norm).sum()
        r2 = rgd.r_g[g2]
        r1 = max(0.6 * r2, waves.rcut)
        g1 = rgd.ceil(r1)
        # Set up confining potential:
        r = rgd.r_g[g1:g2]
        vtr_g = self.vtr_g.copy()
        vtr_g[g1:g2] += scale * np.exp((r2 - r1) / (r1 - r)) / (r - r2)**2
        vtr_g[g2:] = np.inf

        # Nonlocal PAW stuff:
        pt_ng = waves.pt_ng
        dH_nn = waves.dH_nn
        dS_nn = waves.dS_nn
        N = len(pt_ng)

        u_g = rgd.zeros()
        u_ng = rgd.zeros(N)
        duodr_n = np.empty(N)
        a_n = np.empty(N)

        e = waves.e_n[n]
        e0 = e
        ch = Channel(l)
        while True:
            duodr, a = ch.integrate_outwards(u_g, rgd, vtr_g, g1, e)

            for n in range(N):
                duodr_n[n], a_n[n] = ch.integrate_outwards(u_ng[n], rgd,
                                                           vtr_g, g1, e,
                                                           pt_g=pt_ng[n])

            A_nn = (dH_nn - e * dS_nn) / (4 * pi)
            B_nn = rgd.integrate(pt_ng[:, None] * u_ng, -1)
            c_n = rgd.integrate(pt_ng * u_g, -1)
            d_n = np.linalg.solve(np.dot(A_nn, B_nn) + np.eye(N),
                                  np.dot(A_nn, c_n))
            u_g[:g1 + 1] -= np.dot(d_n, u_ng[:, :g1 + 1])
            a -= np.dot(d_n, a_n)
            duodr -= np.dot(duodr_n, d_n)
            uo = u_g[g1]

            duidr = ch.integrate_inwards(u_g, rgd, vtr_g, g1, e, gmax=g2)
            ui = u_g[g1]
            A = duodr / uo - duidr / ui
            u_g[g1:] *= uo / ui
            x = (norm / rgd.integrate(u_g**2, -2) * (4 * pi))**0.5
            u_g *= x
            a *= x

            if abs(A) < 1e-5:
                break

            e += 0.5 * A * u_g[g1]**2

        u_g[1:] /= rgd.r_g[1:]
        u_g[0] = a * 0.0**l
        return u_g, r1, r2, e - e0

    def logarithmic_derivative(self, l, energies, rcut):
        rgd = self.rgd
        ch = Channel(l)
        gcut = rgd.round(rcut)

        N = 0
        if l < len(self.waves_l):
            # Nonlocal PAW stuff:
            waves = self.waves_l[l]
            if len(waves) > 0:
                pt_ng = waves.pt_ng
                dH_nn = waves.dH_nn
                dS_nn = waves.dS_nn
                N = len(pt_ng)

        u_g = rgd.zeros()
        u_ng = rgd.zeros(N)
        dudr_n = np.empty(N)

        logderivs = []
        for e in energies:
            dudr = ch.integrate_outwards(u_g, rgd, self.vtr_g, gcut, e)[0]
            u = u_g[gcut]

            if N:
                for n in range(N):
                    dudr_n[n] = ch.integrate_outwards(u_ng[n], rgd,
                                                      self.vtr_g, gcut, e,
                                                      pt_g=pt_ng[n])[0]

                A_nn = (dH_nn - e * dS_nn) / (4 * pi)
                B_nn = rgd.integrate(pt_ng[:, None] * u_ng, -1)
                c_n = rgd.integrate(pt_ng * u_g, -1)
                d_n = np.linalg.solve(np.dot(A_nn, B_nn) + np.eye(N),
                                      np.dot(A_nn, c_n))
                u -= np.dot(u_ng[:, gcut], d_n)
                dudr -= np.dot(dudr_n, d_n)

            logderivs.append(dudr / u)

        return logderivs

    def make_paw_setup(self, tag=None):
        aea = self.aea

        from gpaw.setup_data import SetupData
        setup = SetupData(aea.symbol, aea.xc.name, tag, readxml=False)

        setup.id_j = []

        J = []  # new reordered j and i indices
        I = []

        # Bound states:
        j = 0
        i = 0
        for l, waves in enumerate(self.waves_l):
            for n, f, e, phi_g, phit_g, pt_g in zip(waves.n_n, waves.f_n,
                                                    waves.e_n, waves.phi_ng,
                                                    waves.phit_ng,
                                                    waves.pt_ng):
                if n != -1:
                    setup.append(n, l, f, e, waves.rcut, phi_g, phit_g, pt_g)
                    id = '%d%s' % (n, 'spdf'[l])
                    setup.id_j.append(id)
                    J.append(j)
                    I.extend(range(i, i + 2 * l + 1))
                j += 1
                i += 2 * l + 1

        # Excited states:
        j = 0
        i = 0
        for l, waves in enumerate(self.waves_l):
            ne = 0
            for n, f, e, phi_g, phit_g, pt_g in zip(waves.n_n, waves.f_n,
                                                    waves.e_n, waves.phi_ng,
                                                    waves.phit_ng,
                                                    waves.pt_ng):
                if n == -1:
                    setup.append(n, l, f, e, waves.rcut, phi_g, phit_g, pt_g)
                    ne += 1
                    id = '%s%d' % ('spdf'[l], ne)
                    setup.id_j.append(id)
                    J.append(j)
                    I.extend(range(i, i + 2 * l + 1))
                j += 1
                i += 2 * l + 1

        nj = sum(len(waves) for waves in self.waves_l)
        e_kin_jj = np.zeros((nj, nj))
        j1 = 0
        for waves in self.waves_l:
            j2 = j1 + len(waves)
            e_kin_jj[j1:j2, j1:j2] = waves.dekin_nn
            j1 = j2
        setup.e_kin_jj = e_kin_jj[J][:, J].copy()

        setup.nc_g = self.nc_g * sqrt(4 * pi)
        setup.nct_g = self.nct_g * sqrt(4 * pi)
        setup.e_kinetic_core = self.ekincore
        setup.vbar_g = self.v0r_g * sqrt(4 * pi)
        setup.vbar_g[1:] /= self.rgd.r_g[1:]
        setup.vbar_g[0] = setup.vbar_g[1]
        setup.Z = aea.Z
        setup.Nc = self.ncore
        setup.Nv = self.nvalence
        setup.e_kinetic = aea.ekin
        setup.e_xc = aea.exc
        setup.e_electrostatic = aea.eH + aea.eZ
        setup.e_total = aea.exc + aea.ekin + aea.eH + aea.eZ
        setup.rgd = self.rgd
        setup.rcgauss = 1 / sqrt(self.alpha)

        self.calculate_exx_integrals()
        setup.ExxC = self.exxcc
        setup.X_p = pack2(self.exxcv_ii[I][:, I])

        setup.tauc_g = self.tauc_g * (4 * pi)**0.5
        setup.tauct_g = self.tauct_g * (4 * pi)**0.5

        if self.aea.scalar_relativistic:
            reltype = 'scalar-relativistic'
        else:
            reltype = 'non-relativistic'
        attrs = [('type', reltype),
                 ('version', 2),
                 ('name', 'gpaw-%s' % version)]
        setup.generatorattrs = attrs

        setup.l0 = self.l0
        setup.e0 = 0.0
        setup.r0 = self.r0
        setup.nderiv0 = self.nderiv0

        setup.basis = self.basis

        if self.core_hole:
            n, l, occ = self.core_hole
            phi_g = self.aea.channels[l].phi_ng[n - l - 1]
            setup.ncorehole = n
            setup.lcorehole = l
            setup.fcorehole = occ
            setup.phicorehole_g = phi_g
            setup.has_corehole = True

        return setup

    def calculate_exx_integrals(self):
        # Find core states:
        core = []
        lmax = 0
        for l, ch in enumerate(self.aea.channels):
            for n, phi_g in enumerate(ch.phi_ng):
                if (l >= len(self.waves_l) or
                    (l < len(self.waves_l) and
                     n + l + 1 not in self.waves_l[l].n_n)):
                    core.append((l, phi_g))
                    if l > lmax:
                        lmax = l

        lmax = max(lmax, len(self.waves_l) - 1)
        G_LLL = make_gaunt(lmax)

        # Calculate core contribution to EXX energy:
        self.exxcc = 0.0
        j1 = 0
        for l1, phi1_g in core:
            f = 1.0
            for l2, phi2_g in core[j1:]:
                n_g = phi1_g * phi2_g
                for l in range((l1 + l2) % 2, l1 + l2 + 1, 2):
                    G = (G_LLL[l1**2:(l1 + 1)**2,
                               l2**2:(l2 + 1)**2,
                               l**2:(l + 1)**2]**2).sum()
                    vr_g = self.rgd.poisson(n_g, l)
                    e = f * self.rgd.integrate(vr_g * n_g, -1) / 4 / pi
                    self.exxcc -= e * G
                f = 2.0
            j1 += 1

        self.log('EXX (core-core):', self.exxcc, 'Hartree')

        # Calculate core-valence contribution to EXX energy:
        ni = sum(len(waves) * (2 * l + 1)
                 for l, waves in enumerate(self.waves_l))

        self.exxcv_ii = np.zeros((ni, ni))

        i1 = 0
        for l1, waves1 in enumerate(self.waves_l):
            for phi1_g in waves1.phi_ng:
                i2 = 0
                for l2, waves2 in enumerate(self.waves_l):
                    for phi2_g in waves2.phi_ng:
                        X_mm = self.exxcv_ii[i1:i1 + 2 * l1 + 1,
                                             i2:i2 + 2 * l2 + 1]
                        if (l1 + l2) % 2 == 0:
                            for lc, phi_g in core:
                                n_g = phi1_g * phi_g
                                for l in range((l1 + lc) % 2,
                                               max(l1, l2) + lc + 1, 2):
                                    vr_g = self.rgd.poisson(phi2_g * phi_g, l)
                                    e = (self.rgd.integrate(vr_g * n_g, -1) /
                                         (4 * pi))
                                    for mc in range(2 * lc + 1):
                                        for m in range(2 * l + 1):
                                            G_L = G_LLL[:,
                                                        lc**2 + mc,
                                                        l**2 + m]
                                            X_mm += np.outer(
                                                G_L[l1**2:(l1 + 1)**2],
                                                G_L[l2**2:(l2 + 1)**2]) * e
                        i2 += 2 * l2 + 1
                i1 += 2 * l1 + 1


def main(argv=None):
    from optparse import OptionParser

    parser = OptionParser(usage='gwap dataset [options] element')
    add = parser.add_option
    add('-f', '--xc-functional', type='string', default='LDA',
        help='Exchange-Correlation functional (default value LDA)',
        metavar='<XC>')
    add('-C', '--configuration',
        help='e.g. for Li: "[(1, 0, 2, -1.878564), (2, 0, 1, -0.10554),'
        '(2, 1, 0, 0.0)]"')
    add('-P', '--projectors',
        help='Projector functions - use comma-separated - ' +
        'nl values, where n can be pricipal quantum number ' +
        '(integer) or energy (floating point number). ' +
        'Example: 2s,0.5s,2p,0.5p,0.0d.')
    add('-r', '--radius',
        help='1.2 or 1.2,1.1,1.1')
    add('-0', '--zero-potential',
        metavar='nderivs,radius',
        help='Parameters for zero potential.')
    add('-c', '--pseudo-core-density-radius', type=float,
        metavar='radius',
        help='Radius for pseudizing core density.')
    add('-z', '--pseudize',
        metavar='type,nderivs',
        help='Parameters for pseudizing wave functions.')
    add('-p', '--plot', action='store_true')
    add('-l', '--logarithmic-derivatives',
        metavar='spdfg,e1:e2:de,radius',
        help='Plot logarithmic derivatives. ' +
        'Example: -l spdf,-1:1:0.05,1.3. ' +
        'Energy range and/or radius can be left out.')
    add('-w', '--write', action='store_true')
    add('-s', '--scalar-relativistic', action='store_true')
    add('--no-check', action='store_true')
    add('-t', '--tag', type='string')
    add('-a', '--alpha', type=float)
    add('-b', '--create-basis-set', action='store_true')
    add('--core-hole')
    add('-e', '--electrons', type=int)
    add('-o', '--output')
    
    opt, symbols = parser.parse_args(argv)

    for symbol in symbols:
        kwargs = get_parameters(symbol, opt)
        gen = _generate(**kwargs)

        if not opt.no_check:
            gen.check_all()

        if opt.create_basis_set or opt.write:
            basis = None  # gen.create_basis_set()
            
            if opt.create_basis_set:
                basis.write_xml()

            if opt.write:
                gen.make_paw_setup(opt.tag).write_xml()

        if opt.logarithmic_derivatives or opt.plot:
            import matplotlib.pyplot as plt
            if opt.logarithmic_derivatives:
                r = 1.1 * gen.rcmax
                emin = min(min(wave.e_n) for wave in gen.waves_l) - 0.8
                emax = max(max(wave.e_n) for wave in gen.waves_l) + 0.8
                lvalues, energies, r = parse_ld_str(
                    opt.logarithmic_derivatives, (emin, emax, 0.05), r)
                ldmax = 0.0
                for l in lvalues:
                    ld = gen.aea.logarithmic_derivative(l, energies, r)
                    plt.plot(energies, ld, colors[l], label='spdfg'[l])
                    ld = gen.logarithmic_derivative(l, energies, r)
                    plt.plot(energies, ld, '--' + colors[l])

                    # Fixed points:
                    if l < len(gen.waves_l):
                        efix = gen.waves_l[l].e_n
                        ldfix = gen.logarithmic_derivative(l, efix, r)
                        plt.plot(efix, ldfix, 'x' + colors[l])
                        ldmax = max(ldmax, max(abs(ld) for ld in ldfix))

                    if l == gen.l0:
                        efix = [0.0]
                        ldfix = gen.logarithmic_derivative(l, efix, r)
                        plt.plot(efix, ldfix, 'x' + colors[l])
                        ldmax = max(ldmax, max(abs(ld) for ld in ldfix))

                if ldmax != 0.0:
                    plt.axis(ymin=-3 * ldmax, ymax=3 * ldmax)
                plt.xlabel('energy [Ha]')
                plt.ylabel(r'$d\phi_{\ell\epsilon}(r)/dr/\phi_{\ell\epsilon}' +
                           r'(r)|_{r=r_c}$')
                plt.legend(loc='best')

            if opt.plot:
                gen.plot()

                if opt.create_basis_set:
                    gen.basis.generatordata = ''  # we already printed this
                    BasisPlotter(show=True).plot(gen.basis)

            try:
                plt.show()
            except KeyboardInterrupt:
                pass

    world.barrier()
    return gen


def get_parameters(symbol, opt):
    if opt.electrons:
        par = parameters[symbol + str(opt.electrons)]
    else:
        Z = atomic_numbers[symbol]
        par = parameters[symbol + str(default[Z])]
        
    projectors, radii = par[:2]
    if len(par) == 3:
        extra = par[2]
    else:
        extra = {}

    if opt.configuration:
        configuration = eval(opt.configuration)
    else:
        configuration = None

    if opt.projectors:
        projectors = opt.projectors

    if opt.radius:
        radii = [float(r) for r in opt.radius.split(',')]

    if isinstance(radii, float):
        radii = [radii]

    if opt.pseudize:
        type, nderiv = opt.pseudize.split(',')
        pseudize = (type, int(nderiv))
    else:
        pseudize = ('poly', 4)

    if opt.zero_potential:
        x = opt.zero_potential.split(',')
        nderiv0 = int(x[0])
        r0 = float(x[1])
    else:
        if projectors[-1].isupper():
            nderiv0 = 5
            r0 = extra.get('r0', min(radii) * 0.9)
        else:
            nderiv0 = 2
            r0 = extra.get('r0', min(radii))

    if opt.pseudo_core_density_radius:
        rcore = opt.pseudo_core_density_radius
    else:
        rcore = extra.get('rcore')
        
    return dict(symbol=symbol,
                xc=opt.xc_functional,
                configuration=configuration,
                projectors=projectors,
                radii=radii,
                scalar_relativistic=opt.scalar_relativistic, alpha=opt.alpha,
                r0=r0, nderiv0=nderiv0,
                pseudize=pseudize, rcore=rcore,
                core_hole=opt.core_hole,
                output=opt.output)


def _generate(symbol, xc, configuration, projectors, radii,
              scalar_relativistic, alpha,
              r0, nderiv0,
              pseudize, rcore, core_hole, output):
    if output is not None:
        output = open(output, 'w')
    aea = AllElectronAtom(symbol, xc, configuration=configuration, log=output)
    gen = PAWSetupGenerator(aea, projectors, scalar_relativistic, core_hole,
                            fd=output)

    gen.construct_shape_function(alpha, radii, eps=1e-10)
    gen.calculate_core_density()
    gen.find_local_potential(r0, nderiv0)
    gen.add_waves(radii)
    gen.pseudize(pseudize[0], pseudize[1], rcore=rcore)
    gen.construct_projectors()
    return gen


def generate_all():
    if len(sys.argv) > 1:
        atoms = sys.argv[1:]
    else:
        atoms = None
    functionals = ['LDA', 'PBE', 'revPBE', 'RPBE']
    for name in parameters:
        n = sum(c.isalpha() for c in name)
        symbol = name[:n]
        electrons = name[n:]
        if atoms and symbol not in atoms:
            continue
        print(name, symbol, electrons)
        for xc in functionals:
            argv = [symbol, '-swf', xc, '-e', electrons, '-t', electrons + 'e']
            if xc == 'PBE':
                argv.append('-b')
            main(argv)


if __name__ == '__main__':
    main()
