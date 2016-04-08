from math import exp, sin, cos, pi, sqrt, acos, asin

import numpy as np

from ase.units import Bohr, Hartree, alpha
import _gpaw
from gpaw.pes import ds_prefactor


class State:

    """Electronic state base class."""

    def __init__(self, gd):
        self.gd = gd
        self.r_c = []
        for c in range(3):
            r_G = (np.arange(gd.n_c[c], dtype=float) +
                   gd.beg_c[c]) * gd.h_cv.diagonal()[c]
            self.r_c.append(r_G)

    def get_grid(self):
        return self.grid_g

    def set_grid(self, grid_g):
#        assert(grid_g.shape == self.gd.empty().shape)
        self.grid_g = grid_g

    def get_energy(self):
        """The states energy in [eV]"""
        return self.energy * Hartree

    def set_energy(self, energy):
        """Set the states energy in [eV]"""
        self.energy = energy / Hartree


class BoundState(State):

    """Bound state originating from a gpaw calculation."""

    def __init__(self,
                 calculator,
                 s, n):
        self.gd = calculator.wfs.gd
        self.kpt = calculator.wfs.kpt_u[s]
        self.grid_g = self.kpt.psit_nG[n]
        self.energy = self.kpt.eps_n[n]


class H1s(State):

    """Analytic H1s state."""

    def __init__(self, gd=None, center=None, Z=1):
        if center is None:
            center = 0.5 * gd.n_c * gd.h_c
        else:
            self.center = center / Bohr
        self.Z = Z
        self.energy = - Z ** 2 / 2.
        if gd is not None:
            self.grid_g = None
            State.__init__(self, gd)

    def get_grid(self):
        if self.grid_g is None:
            gd = self.gd
            assert(gd.orthogonal)
            wf = gd.zeros(dtype=float)
            vr0 = np.array([0., 0., 0.])
            for i in range(gd.n_c[0]):
                vr0[0] = (gd.beg_c[0] + i) * gd.h_cv[0, 0]
                for j in range(gd.n_c[1]):
                    vr0[1] = (gd.beg_c[1] + j) * gd.h_cv[1, 1]
                    for k in range(gd.n_c[2]):
                        vr0[2] = (gd.beg_c[2] + k) * gd.h_cv[2, 2]
                        vr = vr0 - self.center
                        r = sqrt(np.dot(vr, vr))
                        wf[i, j, k] = exp(-self.Z * r)
            self.grid_g = wf * sqrt(self.Z ** 3 / pi)
        return self.grid_g

    def get_me_c(self, k_c, form='L'):
        """Transition matrix element."""
        k = sqrt(np.dot(k_c, k_c))
        if form == 'L':
            pre = 2j
        elif form == 'V':
            pre = 1j
        else:
            raise NonImplementedError
        return pre * k_c * self.FT(k)

    def get_ds(self, Ekin, form='L', units='Mb'):
        """Angular averaged cross section.

        Ekin: photoelectron kinetic energy [eV]"""
        E = Ekin / Hartree
        k = sqrt(2 * E)
        omega = E - self.energy
        k_c = np.array([0., 0., k])
        me_c = self.get_me_c(k_c, form)
        T2 = np.abs(np.dot(me_c, me_c)) / 3.
        pre = ds_prefactor[units]
        return pre * (2 * pi) ** 2 * alpha * k * 4 * pi / omega * T2

    def FT(self, k):
        return sqrt(8 * self.Z ** 5) / pi / (k ** 2 + self.Z ** 2) ** 2
