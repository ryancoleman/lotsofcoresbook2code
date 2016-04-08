"""Appelbaum-Hamann local pseudo-potential for silicon.

See::

  Self-Consistent Pseudopotential for Si
  Joel A. Appelbaum and D. R. Hamann
  PRB 8, 1777 (1973)

"""

import numpy as np

from gpaw.setup import BaseSetup
from gpaw.spline import Spline
from gpaw.basis_data import Basis


class AppelbaumHamann(BaseSetup):
    def __init__(self, alpha=0.6102, v1=3.042, v2=-1.372):
        self.alpha = alpha
        self.v1 = v1
        self.v2 = v2

        self.E = 0.0
        self.Z = 14
        self.Nc = 10
        self.Nv = 4
        self.nao = None
        nullspline = Spline(0, 0.5, [0., 0., 0.])
        self.pt_j = [nullspline]
        self.ni = 1
        self.l_j = [0]
        self.f_j = [4]
        self.n_j = [1]
        self.nct = nullspline
        self.Nct = 0.0
        rc = 4.0
        r2_g = np.linspace(0, rc, 100)**2
        x_g = np.exp(-alpha * r2_g)
        self.ghat_l = [Spline(0, rc, 4 * alpha**1.5 / np.pi**0.5 * x_g)]
        self.vbar = Spline(0, rc, 2 * np.pi**0.5 * (v1 + v2 * r2_g) * x_g)
        self.Delta_pL = np.zeros((1, 1))
        self.Delta0 = -4 / (4 * np.pi)**0.5
        self.lmax = 0
        self.K_p = self.M_p = self.MB_p = np.zeros(1)
        self.M_pp = np.zeros((1, 1))
        self.Kc = 0.0
        self.MB = 0.0
        self.M = 0.0
        self.xc_correction = None
        self.HubU = None
        self.dO_ii = np.zeros((1, 1))
        self.type = 'ah'
        self.fingerprint = None

    def build(self, basis):
        if basis is None:
            basis = Basis('Si', 'sz(dzp)')
        elif isinstance(basis, str):
            basis = Basis('Si', basis)
        self.basis = basis
        self.phit_j = self.basis.tosplines()
        self.nao = self.basis.nao

    def print_info(self, text):
        text('Appelbaum-Hamann pseudo potential')

    def calculate_initial_occupation_numbers(self, magmom, hund, charge,
                                             nspins):
        assert nspins == 1
        return np.array([(2.0, 2.0 / 3, 2.0 / 3, 2.0 / 3)])

    def initialize_density_matrix(self, f_si):
        return np.zeros((len(f_si), 1))

    def get_default_nbands(self):
        return 3
