"""BSSE: Basis Set Superposition Error module.

Defines a Setup-like class which has no properties that change anything,
except for an atomic basis set."""

import numpy as np
from ase.data import atomic_numbers

from gpaw.setup import BaseSetup
from gpaw.setup_data import SetupData
from gpaw.basis_data import Basis
from gpaw.spline import Spline
from gpaw.utilities import min_locfun_radius


# Some splines are mandatory, but should then be zero to avoid affecting things
zero_function = Spline(0, min_locfun_radius, [0.0, 0.0, 0.0])

# Some operations fail horribly if the splines are zero, due to weird
# divisions and assumptions that various quantities are nonzero
#
# We'll use a function which is almost zero for these things
nonzero_function = Spline(0, min_locfun_radius, [0.0, 1.0e-12, 0.0]) # XXX

class GhostSetup(BaseSetup):
    def __init__(self, basis, data):
        self.symbol = data.symbol
        self.data = data
        self.phit_j = basis.tosplines()
        self.basis = basis
        self.nao = sum([2 * phit.get_angular_momentum_number() + 1
                        for phit in self.phit_j])
        self.HubU = None
        self.filename = None
        self.fingerprint = None
        self.type = 'ghost'

        self.Z = 0
        self.Nv = 0
        self.Nc = 0

        self.ni = 1
        self.pt_j = [zero_function]
        self.wg_lg = None
        self.g_lg = None

        self.Nct = 1e-12 # XXX XXX XXX XXX
        self.nct = nonzero_function # XXXXXX
        self.lmax = 0
        self.xc_correction = None
        self.ghat_l = [nonzero_function] * (self.lmax + 1) # XXXXXX
        self.rcgauss = 1e12 # XXX XXX XXX XXX
        self.vbar = zero_function

        self.Delta0 = 0.0
        self.Delta_pL = np.zeros((1, self.lmax + 1))

        self.E = 0.0
        self.Kc = 0.0
        self.M = 0.0
        self.M_p = np.zeros(1)
        self.M_pp = np.zeros((1, 1))
        self.K_p = np.zeros(1)
        self.MB = 0.0
        self.MB_p = np.zeros(1)
        self.dO_ii = np.zeros((1, 1))
        self.f_j = [0.0]
        self.n_j = [0]
        self.l_j = [0]
        self.l_orb_j = [0]
        self.nj = 1
        self.lq = None # XXXX

        self.rcutfilter = None
        self.rcore = None
        self.N0_p = np.zeros(1)
        self.nabla_iiv = None
        self.rnabla_iiv = None
        self.rxp_iiv = None
        self.phicorehole_g = None
        self.rgd = None
        self.rcut_j = [0.5]
        self.tauct = None
        self.Delta_iiL = None
        self.B_ii = None
        self.dC_ii = None
        self.X_p = None
        self.ExxC = None
        self.dEH0 = 0.0
        self.dEH_p = np.zeros(1)
        self.extra_xc_data = {}


class GhostSetupData:
    def __init__(self, symbol):
        self.chemsymbol = symbol
        self.symbol = symbol + '.ghost'
        self.Z = atomic_numbers[symbol]

    def build(self, xcfunc, lmax, basis, filter=None):
        if basis is None:
            raise ValueError('Loading partial waves not supported right now')
        setup = GhostSetup(basis, self)
        return setup

    def print_info(self, text, _setup):
        text('Ghost setup for %s' % self.chemsymbol)
