import numpy as np

from gpaw.atom.atompaw import AtomPAW
from gpaw.utilities import erf
from gpaw.setup import BaseSetup
from gpaw.spline import Spline
from gpaw.basis_data import Basis, BasisFunction

null_spline = Spline(0, 1.0, [0., 0., 0.])


def screen_potential(r, v, charge, rcut=None, a=None):
    """Split long-range potential into short-ranged contributions.

    The potential v is a long-ranted potential with the asymptotic form Z/r
    corresponding to the given charge.
    
    Return a potential vscreened and charge distribution rhocomp such that

      v(r) = vscreened(r) + vHartree[rhocomp](r).

    The returned quantities are truncated to a reasonable cutoff radius.
    """
    vr = v * r + charge
    
    if rcut is None:
        err = 0.0
        i = len(vr)
        while err < 1e-5:
            # Things can be a bit sensitive to the threshold.  The O.pz-mt
            # setup gets 20-30 Bohr long compensation charges if it's 1e-6.
            i -= 1
            err = abs(vr[i])
        i += 1

        icut = np.searchsorted(r, r[i] * 1.1)
    else:
        icut = np.searchsorted(r, rcut)
    rcut = r[icut]
    rshort = r[:icut]

    if a is None:
        a = rcut / 5.0 # XXX why is this so important?
    vcomp = charge * erf(rshort / (np.sqrt(2.0) * a)) / rshort
    # XXX divide by r
    rhocomp = charge * (np.sqrt(2.0 * np.pi) * a)**(-3) * \
        np.exp(-0.5 * (rshort / a)**2)
    vscreened = v[:icut] + vcomp
    return vscreened, rhocomp


def figure_out_valence_states(ppdata):
    from gpaw.atom.configurations import configurations
    from ase.data import chemical_symbols
    # ppdata.symbol may not be a chemical symbol so use Z
    chemical_symbol = chemical_symbols[ppdata.Z]
    Z, config = configurations[chemical_symbol]
    assert Z == ppdata.Z
    
    # Okay, we need to figure out occupations f_ln when we don't know
    # any info about existing states on the pseudopotential.
    #
    # The plan is to loop over all states and count until only the correct
    # number of valence electrons "remain".
    nelectrons = 0
    ncore = ppdata.Z - ppdata.Nv
    assert ppdata.Nv > 0
    iterconfig = iter(config)
    for n, l, occ, eps in iterconfig:
        nelectrons += occ
        if nelectrons == ncore:
            break
        elif nelectrons >= ncore:
            raise ValueError('Cannot figure out what states should exist '
                             'on this pseudopotential.')
    
    f_ln = {}
    l_j = []
    f_j = []
    n_j = []
    for n, l, occ, eps in iterconfig:
        f_ln.setdefault(l, []).append(occ)
        l_j.append(l)
        f_j.append(occ)
        n_j.append(n)
    lmax = max(f_ln.keys())
    f_ln = [f_ln.get(l, []) for l in range(lmax + 1)]
    return n_j, l_j, f_j, f_ln


def generate_basis_functions(ppdata):
    class SimpleBasis(Basis):
        def __init__(self, symbol, l_j):
            Basis.__init__(self, symbol, 'simple', readxml=False)
            self.generatordata = 'simple'
            self.d = 0.02
            self.ng = 160
            rgd = self.get_grid_descriptor()
            bf_j = self.bf_j
            rcgauss = rgd.r_g[-1] / 3.0
            gauss_g = np.exp(-(rgd.r_g / rcgauss)**2.0)
            for l in l_j:
                phit_g = rgd.r_g**l * gauss_g
                norm = (rgd.integrate(phit_g**2) / (4 * np.pi))**0.5
                phit_g /= norm
                bf = BasisFunction(l, rgd.r_g[-1], phit_g, 'gaussian')
                bf_j.append(bf)
    #l_orb_j = [state.l for state in self.data['states']]
    b1 = SimpleBasis(ppdata.symbol, ppdata.l_orb_j)
    apaw = AtomPAW(ppdata.symbol, [ppdata.f_ln], h=0.05, rcut=9.0,
                   basis={ppdata.symbol: b1},
                   setups={ppdata.symbol: ppdata},
                   lmax=0, txt=None)
    basis = apaw.extract_basis_functions()
    return basis


def pseudoplot(pp, show=True):
    import pylab as pl
    
    fig = pl.figure()
    wfsax = fig.add_subplot(221)
    ptax = fig.add_subplot(222)
    vax = fig.add_subplot(223)
    rhoax = fig.add_subplot(224)

    def spline2grid(spline):
        rcut = spline.get_cutoff()
        r = np.linspace(0.0, rcut, 2000)
        return r, spline.map(r)

    for phit in pp.phit_j:
        r, y = spline2grid(phit)
        wfsax.plot(r, y, label='wf l=%d' % phit.get_angular_momentum_number())

    for pt in pp.pt_j:
        r, y = spline2grid(pt)
        ptax.plot(r, y, label='pr l=%d' % pt.get_angular_momentum_number())

    for ghat in pp.ghat_l:
        r, y = spline2grid(ghat)
        rhoax.plot(r, y, label='cc l=%d' % ghat.get_angular_momentum_number())

    r, y = spline2grid(pp.vbar)
    vax.plot(r, y, label='vbar')
    
    vax.set_ylabel('potential')
    rhoax.set_ylabel('density')
    wfsax.set_ylabel('wfs')
    ptax.set_ylabel('projectors')

    for ax in [vax, rhoax, wfsax, ptax]:
        ax.legend()

    if show:
        pl.show()


class PseudoPotential(BaseSetup):
    def __init__(self, data, basis=None):
        self.data = data

        self.R_sii = None
        self.HubU = None
        self.lq = None

        self.filename = None
        self.fingerprint = None
        self.symbol = data.symbol
        self.type = data.name

        self.Z = data.Z
        self.Nv = data.Nv
        self.Nc = data.Nc

        self.f_j = data.f_j
        self.n_j = data.n_j
        self.l_j = data.l_j
        self.l_orb_j = data.l_orb_j
        self.nj = len(data.l_j)

        self.ni = sum([2 * l + 1 for l in data.l_j])
        self.pt_j = data.get_projectors()
        if len(self.pt_j) == 0:
            assert False # not sure yet about the consequences of
            # cleaning this up in the other classes
            self.l_j = [0]
            self.pt_j = [null_spline]
        
        if basis is None:
            basis = data.create_basis_functions()
        self.phit_j = basis.tosplines()
        self.basis = basis
        self.nao = sum([2 * phit.get_angular_momentum_number() + 1
                        for phit in self.phit_j])

        self.Nct = 0.0
        self.nct = null_spline

        self.lmax = 0

        self.xc_correction = None

        r, l_comp, g_comp = data.get_compensation_charge_functions()
        self.ghat_l = [Spline(l, r[-1], g) for l, g in zip(l_comp, g_comp)]
        self.rcgauss = data.rcgauss

        # accuracy is rather sensitive to this
        self.vbar = data.get_local_potential()

        _np = self.ni * (self.ni + 1) // 2
        self.Delta0 = data.Delta0
        self.Delta_pL = np.zeros((_np, 1))

        self.E = 0.0
        self.Kc = 0.0
        self.M = 0.0
        self.M_p = np.zeros(_np)
        self.M_pp = np.zeros((_np, _np))

        self.K_p = data.expand_hamiltonian_matrix()
        self.MB = 0.0
        self.MB_p = np.zeros(_np)
        self.dO_ii = np.zeros((self.ni, self.ni))

        # We don't really care about these variables
        self.rcutfilter = None
        self.rcore = None

        self.N0_p = np.zeros(_np) # not really implemented
        self.nabla_iiv = None
        self.rnabla_iiv = None
        self.rxp_iiv = None
        self.phicorehole_g = None
        self.rgd = data.rgd
        self.rcut_j = data.rcut_j
        self.tauct = None
        self.Delta_iiL = None
        self.B_ii = None
        self.dC_ii = None
        self.X_p = None
        self.ExxC = None
        self.dEH0 = 0.0
        self.dEH_p = np.zeros(_np)
        self.extra_xc_data = {}

        self.wg_lg = None
        self.g_lg = None
