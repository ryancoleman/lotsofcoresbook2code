# Todo: ACDF formula
from __future__ import division

import sys
from math import pi

import numpy as np
from ase.units import Hartree
from ase.utils import prnt

import gpaw.mpi as mpi
from gpaw.xc import XC
from gpaw.xc.tools import vxc
from gpaw.xc.kernel import XCNull
from gpaw.response.pair import PairDensity
from gpaw.wavefunctions.pw import PWDescriptor
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.utilities import unpack, unpack2, packed_index
from gpaw.response.wstc import WignerSeitzTruncatedCoulomb


def pawexxvv(atomdata, D_ii):
    """PAW correction for valence-valence EXX energy."""
    ni = len(D_ii)
    V_ii = np.empty((ni, ni))
    for i1 in range(ni):
        for i2 in range(ni):
            V = 0.0
            for i3 in range(ni):
                p13 = packed_index(i1, i3, ni)
                for i4 in range(ni):
                    p24 = packed_index(i2, i4, ni)
                    V += atomdata.M_pp[p13, p24] * D_ii[i3, i4]
            V_ii[i1, i2] = V
    return V_ii
    
    
def select_kpts(kpts, calc):
    """Function to process input parameters that take a list of k-points given
    in different format and returns a list of indices of the corresponding
    k-points in the IBZ."""
    if kpts is None:
        # Do all k-points in the IBZ:
        return range(calc.wfs.kd.nibzkpts)
    
    if np.asarray(kpts).ndim == 1:
        return kpts
    
    # Find k-points:
    bzk_Kc = calc.get_bz_k_points()
    indices = []
    for k_c in kpts:
        d_Kc = bzk_Kc - k_c
        d_Kc -= d_Kc.round()
        K = abs(d_Kc).sum(1).argmin()
        if not np.allclose(d_Kc[K], 0):
            raise ValueError('Could not find k-point: {0}'.format(k_c))
        k = calc.wfs.kd.bz2ibz_k[K]
        indices.append(k)
    return indices
    
    
class EXX(PairDensity):
    def __init__(self, calc, xc=None, kpts=None, bands=None, ecut=None,
                 omega=None, world=mpi.world, txt=sys.stdout, timer=None):
    
        PairDensity.__init__(self, calc, ecut, world=world, txt=txt,
                             timer=timer)

        if xc is None or xc == 'EXX':
            self.exx_fraction = 1.0
            xc = XC(XCNull())
        elif xc == 'PBE0':
            self.exx_fraction = 0.25
            xc = XC('HYB_GGA_XC_PBEH')
        elif xc == 'HSE03':
            omega = 0.106
            self.exx_fraction = 0.25
            xc = XC('HYB_GGA_XC_HSE03')
        elif xc == 'HSE06':
            omega = 0.11
            self.exx_fraction = 0.25
            xc = XC('HYB_GGA_XC_HSE06')
        elif xc == 'B3LYP':
            self.exx_fraction = 0.2
            xc = XC('HYB_GGA_XC_B3LYP')
            
        self.xc = xc
        self.omega = omega
        self.exc = np.nan  # density dependent part of xc-energy
        
        self.kpts = select_kpts(kpts, self.calc)
        
        if bands is None:
            # Do all occupied bands:
            bands = [0, self.nocc2]
        
        prnt('Calculating exact exchange contributions for band index',
             '%d-%d' % (bands[0], bands[1] - 1), file=self.fd)
        prnt('for IBZ k-points with indices:',
             ', '.join(str(i) for i in self.kpts), file=self.fd)
        
        self.bands = bands

        if self.ecut is None:
            self.ecut = self.calc.wfs.pd.ecut
        prnt('Plane-wave cutoff: %.3f eV' % (self.ecut * Hartree),
             file=self.fd)
        
        shape = (self.calc.wfs.nspins, len(self.kpts), bands[1] - bands[0])
        self.exxvv_sin = np.zeros(shape)   # valence-valence exchange energies
        self.exxvc_sin = np.zeros(shape)   # valence-core exchange energies
        self.f_sin = np.empty(shape)       # occupation numbers

        # The total EXX energy will not be calculated if we are only
        # interested in a few eigenvalues for a few k-points
        self.exx = np.nan    # total EXX energy
        self.exxvv = np.nan  # valence-valence
        self.exxvc = np.nan  # valence-core
        self.exxcc = 0.0     # core-core

        self.mysKn1n2 = None  # my (s, K, n1, n2) indices
        self.distribute_k_points_and_bands(0, self.nocc2)
        
        # All occupied states:
        self.mykpts = [self.get_k_point(s, K, n1, n2)
                       for s, K, n1, n2 in self.mysKn1n2]

        prnt('Using Wigner-Seitz truncated coulomb interaction.',
             file=self.fd)
        self.wstc = WignerSeitzTruncatedCoulomb(self.calc.wfs.gd.cell_cv,
                                                self.calc.wfs.kd.N_c,
                                                self.fd)
        self.iG_qG = {}  # cache
            
        # PAW matrices:
        self.V_asii = []  # valence-valence correction
        self.C_aii = []   # valence-core correction
        self.initialize_paw_exx_corrections()
        
    def calculate(self):
        kd = self.calc.wfs.kd
        nspins = self.calc.wfs.nspins
        
        for s in range(nspins):
            for i, k1 in enumerate(self.kpts):
                K1 = kd.ibz2bz_k[k1]
                kpt1 = self.get_k_point(s, K1, *self.bands)
                self.f_sin[s, i] = kpt1.f_n
                for kpt2 in self.mykpts:
                    if kpt2.s == s:
                        self.calculate_q(i, kpt1, kpt2)
                
                self.calculate_paw_exx_corrections(i, kpt1)

        self.world.sum(self.exxvv_sin)
        
        # Calculate total energy if we have everything needed:
        if (len(self.kpts) == kd.nibzkpts and
            self.bands[0] == 0 and
            self.bands[1] >= self.nocc2):
            exxvv_i = (self.exxvv_sin * self.f_sin).sum(axis=2).sum(axis=0)
            exxvc_i = 2 * (self.exxvc_sin * self.f_sin).sum(axis=2).sum(axis=0)
            self.exxvv = np.dot(kd.weight_k[self.kpts], exxvv_i) / nspins
            self.exxvc = np.dot(kd.weight_k[self.kpts], exxvc_i) / nspins
            self.exx = self.exxvv + self.exxvc + self.exxcc
            prnt('Exact exchange energy:', file=self.fd)
            for kind, exx in [('valence-valence', self.exxvv),
                              ('valence-core', self.exxvc),
                              ('core-core', self.exxcc),
                              ('total', self.exx)]:
                prnt('%16s%11.3f eV' % (kind + ':', exx * Hartree),
                     file=self.fd)
            
            self.exc = self.calculate_hybrid_correction()

        exx_sin = self.exxvv_sin + self.exxvc_sin
        prnt('EXX eigenvalue contributions in eV:', file=self.fd)
        prnt(np.array_str(exx_sin * Hartree, precision=3), file=self.fd)
    
    def get_exx_energy(self):
        return self.exx * Hartree
    
    def get_total_energy(self):
        ham = self.calc.hamiltonian
        return (self.exx * self.exx_fraction + self.exc +
                ham.Etot - ham.Exc) * Hartree
        
    def get_eigenvalue_contributions(self):
        if self.reader is not None:
            self.calc.wfs.read_projections(self.reader)
        b1, b2 = self.bands
        e_sin = vxc(self.calc, self.xc)[:, self.kpts, b1:b2] / Hartree
        e_sin += (self.exxvv_sin + self.exxvc_sin) * self.exx_fraction
        return e_sin * Hartree
        
    def calculate_q(self, i, kpt1, kpt2):
        wfs = self.calc.wfs
        
        q_c = wfs.kd.bzk_kc[kpt2.K] - wfs.kd.bzk_kc[kpt1.K]
        qd = KPointDescriptor([q_c])
        pd = PWDescriptor(self.ecut, wfs.gd, wfs.dtype, kd=qd)

        Q_G = self.get_fft_indices(kpt1.K, kpt2.K, q_c, pd,
                                   kpt1.shift_c - kpt2.shift_c)

        Q_aGii = self.initialize_paw_corrections(pd, soft=True)
        
        for n in range(kpt1.n2 - kpt1.n1):
            ut1cc_R = kpt1.ut_nR[n].conj()
            C1_aGi = [np.dot(Q_Gii, P1_ni[n].conj())
                      for Q_Gii, P1_ni in zip(Q_aGii, kpt1.P_ani)]
            n_mG = self.calculate_pair_densities(ut1cc_R, C1_aGi, kpt2,
                                                 pd, Q_G)
            e = self.calculate_n(pd, n, n_mG, kpt2)
            self.exxvv_sin[kpt1.s, i, n] += e

    def calculate_n(self, pd, n, n_mG, kpt2):
        iG_G = self.get_coulomb_kernel(pd)
        x = 4 * pi / self.calc.wfs.kd.nbzkpts / pd.gd.dv**2

        e = 0.0
        for f, n_G in zip(kpt2.f_n, n_mG):
            x_G = n_G * iG_G
            e -= x * f * pd.integrate(x_G, x_G).real

        return e

    def get_coulomb_kernel(self, pd):
        if self.omega is not None:
            G2_G = pd.G2_qG[0]
            iG_G = np.empty_like(G2_G)
            if pd.kd.gamma:
                iG_G[0] = 1 / (2 * self.omega)
            else:
                iG_G[0] = ((1 - np.exp(-G2_G[0] / (4 * self.omega**2))) /
                           G2_G[0])**0.5
            iG_G[1:] = ((1 - np.exp(-G2_G[1:] / (4 * self.omega**2))) /
                        G2_G[1:])**0.5
            return iG_G

        key = tuple((pd.kd.bzk_kc[0] * self.calc.wfs.kd.N_c).round())
        iG_G = self.iG_qG.get(key)
        if iG_G is None:
            v_G = self.wstc.get_potential(pd)
            iG_G = (v_G / (4 * pi))**0.5
            self.iG_qG[key] = iG_G
        return iG_G

    def initialize_paw_exx_corrections(self):
        for a, atomdata in enumerate(self.calc.wfs.setups):
            V_sii = []
            for D_p in self.calc.density.D_asp[a]:
                D_ii = unpack2(D_p)
                V_ii = pawexxvv(atomdata, D_ii)
                V_sii.append(V_ii)
            C_ii = unpack(atomdata.X_p)
            self.V_asii.append(V_sii)
            self.C_aii.append(C_ii)
            self.exxcc += atomdata.ExxC

    def calculate_paw_exx_corrections(self, i, kpt):
        x = self.calc.wfs.nspins / self.world.size
        s = kpt.s
        
        for V_sii, C_ii, P_ni in zip(self.V_asii, self.C_aii, kpt.P_ani):
            V_ii = V_sii[s]
            v_n = (np.dot(P_ni, V_ii) * P_ni.conj()).sum(axis=1).real
            c_n = (np.dot(P_ni, C_ii) * P_ni.conj()).sum(axis=1).real
            self.exxvv_sin[s, i] -= v_n * x
            self.exxvc_sin[s, i] -= c_n

    def calculate_hybrid_correction(self):
        dens = self.calc.density
        if dens.nt_sg is None:
            dens.interpolate_pseudo_density()
        exc = self.xc.calculate(dens.finegd, dens.nt_sg)
        for a, D_sp in dens.D_asp.items():
            atomdata = dens.setups[a]
            exc += self.xc.calculate_paw_correction(atomdata, D_sp)
        return exc
