# Copyright (C) 2010  CAMd
# Please see the accompanying LICENSE file for further information.

"""This module provides all the classes and functions associated with the
evaluation of exact exchange with k-point sampling."""
from __future__ import print_function
from time import time
from math import pi, sqrt

import numpy as np
from ase.utils import prnt
from ase.units import Hartree
from ase.dft.kpoints import monkhorst_pack
from ase.utils.timing import Timer

import gpaw.mpi as mpi
import gpaw.fftw as fftw
from gpaw.xc.hybrid import HybridXCBase
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.wavefunctions.pw import PWDescriptor, PWLFC
from gpaw.utilities import pack, unpack2, packed_index, logfile, erf
from gpaw.utilities.ewald import madelung


class HybridXC(HybridXCBase):
    orbital_dependent = True

    def __init__(self, name, hybrid=None, xc=None,
                 alpha=None,
                 gamma_point=1,
                 method='standard',
                 bandstructure=False,
                 logfilename='-', bands=None,
                 fcut=1e-10,
                 molecule=False,
                 qstride=1,
                 world=None):
        """Mix standard functionals with exact exchange.

        name: str
            Name of functional: EXX, PBE0, HSE03, HSE06
        hybrid: float
            Fraction of exact exchange.
        xc: str or XCFunctional object
            Standard DFT functional with scaled down exchange.
        method: str
            Use 'standard' standard formula and 'acdf for
            adiabatic-connection dissipation fluctuation formula.
        alpha: float
            XXX describe
        gamma_point: bool
            0: Skip k2-k1=0 interactions.
            1: Use the alpha method.
            2: Integrate the gamma point.
        bandstructure: bool
            Calculate bandstructure instead of just the total energy.
        bands: list of int
            List of bands to calculate bandstructure for.  Default is
            all bands.
        molecule: bool
            Decouple electrostatic interactions between periodically
            repeated images.
        fcut: float
            Threshold for empty band.
        """

        self.alpha = alpha
        self.fcut = fcut

        self.gamma_point = gamma_point
        self.method = method
        self.bandstructure = bandstructure
        self.bands = bands

        self.fd = logfilename
        self.write_timing_information = True

        HybridXCBase.__init__(self, name, hybrid, xc)

        # EXX energies:
        self.exx = None  # total
        self.evv = None  # valence-valence (pseudo part)
        self.evvacdf = None  # valence-valence (pseudo part)
        self.devv = None  # valence-valence (PAW correction)
        self.evc = None  # valence-core
        self.ecc = None  # core-core

        self.exx_skn = None  # bandstructure

        self.qlatest = None

        if world is None:
            world = mpi.world
        self.world = world

        self.molecule = molecule
        
        if isinstance(qstride, int):
            qstride = [qstride] * 3
        self.qstride_c = np.asarray(qstride)
        
        self.timer = Timer()

    def log(self, *args, **kwargs):
        prnt(file=self.fd, *args, **kwargs)
        self.fd.flush()

    def calculate_radial(self, rgd, n_sLg, Y_L, v_sg,
                         dndr_sLg=None, rnablaY_Lv=None,
                         tau_sg=None, dedtau_sg=None):
        return self.xc.calculate_radial(rgd, n_sLg, Y_L, v_sg,
                                        dndr_sLg, rnablaY_Lv)
    
    def calculate_paw_correction(self, setup, D_sp, dEdD_sp=None,
                                 addcoredensity=True, a=None):
        return self.xc.calculate_paw_correction(setup, D_sp, dEdD_sp,
                                 addcoredensity, a)
    
    def initialize(self, dens, ham, wfs, occupations):
        assert wfs.bd.comm.size == 1

        self.xc.initialize(dens, ham, wfs, occupations)

        self.dens = dens
        self.wfs = wfs

        # Make a k-point descriptor that is not distributed
        # (self.kd.comm is serial_comm):
        self.kd = wfs.kd.copy()

        self.fd = logfile(self.fd, self.world.rank)

        wfs.initialize_wave_functions_from_restart_file()

    def set_positions(self, spos_ac):
        self.spos_ac = spos_ac

    def calculate(self, gd, n_sg, v_sg=None, e_g=None):
        # Normal XC contribution:
        exc = self.xc.calculate(gd, n_sg, v_sg, e_g)

        # Add EXX contribution:
        return exc + self.exx * self.hybrid

    def calculate_exx(self):
        """Non-selfconsistent calculation."""

        self.timer.start('EXX')
        self.timer.start('Initialization')
        
        kd = self.kd
        wfs = self.wfs

        if fftw.FFTPlan is fftw.NumpyFFTPlan:
            self.log('NOT USING FFTW !!')

        self.log('Spins:', self.wfs.nspins)

        W = max(1, self.wfs.kd.comm.size // self.wfs.nspins)
        # Are the k-points distributed?
        kparallel = (W > 1)

        # Find number of occupied bands:
        self.nocc_sk = np.zeros((self.wfs.nspins, kd.nibzkpts), int)
        for kpt in self.wfs.kpt_u:
            for n, f in enumerate(kpt.f_n):
                if abs(f) < self.fcut:
                    self.nocc_sk[kpt.s, kpt.k] = n
                    break
            else:
                self.nocc_sk[kpt.s, kpt.k] = self.wfs.bd.nbands
        self.wfs.kd.comm.sum(self.nocc_sk)

        noccmin = self.nocc_sk.min()
        noccmax = self.nocc_sk.max()
        self.log('Number of occupied bands (min, max): %d, %d' %
                 (noccmin, noccmax))
        
        self.log('Number of valence electrons:', self.wfs.setups.nvalence)

        if self.bandstructure:
            self.log('Calculating eigenvalue shifts.')

            # allocate array for eigenvalue shifts:
            self.exx_skn = np.zeros((self.wfs.nspins,
                                     kd.nibzkpts,
                                     self.wfs.bd.nbands))

            if self.bands is None:
                noccmax = self.wfs.bd.nbands
            else:
                noccmax = max(max(self.bands) + 1, noccmax)

        N_c = self.kd.N_c

        vol = wfs.gd.dv * wfs.gd.N_c.prod()
        if self.alpha is None:
            alpha = 6 * vol**(2 / 3.0) / pi**2
        else:
            alpha = self.alpha
        if self.gamma_point == 1:
            if alpha == 0.0:
                qvol = (2*np.pi)**3 / vol / N_c.prod()
                self.gamma = 4*np.pi * (3*qvol / (4*np.pi))**(1/3.) / qvol
            else:
                self.gamma = self.calculate_gamma(vol, alpha)
        else:
            kcell_cv = wfs.gd.cell_cv.copy()
            kcell_cv[0] *= N_c[0]
            kcell_cv[1] *= N_c[1]
            kcell_cv[2] *= N_c[2]
            self.gamma = madelung(kcell_cv) * vol * N_c.prod() / (4 * np.pi)

        self.log('Value of alpha parameter: %.3f Bohr^2' % alpha)
        self.log('Value of gamma parameter: %.3f Bohr^2' % self.gamma)
            
        # Construct all possible q=k2-k1 vectors:
        Nq_c = (N_c - 1) // self.qstride_c
        i_qc = np.indices(Nq_c * 2 + 1, float).transpose(
            (1, 2, 3, 0)).reshape((-1, 3))
        self.bzq_qc = (i_qc - Nq_c) / N_c * self.qstride_c
        self.q0 = ((Nq_c * 2 + 1).prod() - 1) // 2  # index of q=(0,0,0)
        assert not self.bzq_qc[self.q0].any()

        # Count number of pairs for each q-vector:
        self.npairs_q = np.zeros(len(self.bzq_qc), int)
        for s in range(kd.nspins):
            for k1 in range(kd.nibzkpts):
                for k2 in range(kd.nibzkpts):
                    for K2, q, n1_n, n2 in self.indices(s, k1, k2):
                        self.npairs_q[q] += len(n1_n)

        self.npairs0 = self.npairs_q.sum()  # total number of pairs

        self.log('Number of pairs:', self.npairs0)

        # Distribute q-vectors to Q processors:
        Q = self.world.size // self.wfs.kd.comm.size
        myrank = self.world.rank // self.wfs.kd.comm.size
        rank = 0
        N = 0
        myq = []
        nq = 0
        for q, n in enumerate(self.npairs_q):
            if n > 0:
                nq += 1
                if rank == myrank:
                    myq.append(q)
            N += n
            if N >= (rank + 1.0) * self.npairs0 / Q:
                rank += 1

        assert len(myq) > 0, 'Too few q-vectors for too many processes!'
        self.bzq_qc = self.bzq_qc[myq]
        try:
            self.q0 = myq.index(self.q0)
        except ValueError:
            self.q0 = None

        self.log('%d x %d x %d k-points' % tuple(self.kd.N_c))
        self.log('Distributing %d IBZ k-points over %d process(es).' %
                 (kd.nibzkpts, self.wfs.kd.comm.size))
        self.log('Distributing %d q-vectors over %d process(es).' % (nq, Q))

        # q-point descriptor for my q-vectors:
        qd = KPointDescriptor(self.bzq_qc)

        # Plane-wave descriptor for all wave-functions:
        self.pd = PWDescriptor(wfs.pd.ecut, wfs.gd,
                               dtype=wfs.pd.dtype, kd=kd)

        # Plane-wave descriptor pair-densities:
        self.pd2 = PWDescriptor(self.dens.pd2.ecut, self.dens.gd,
                                dtype=wfs.dtype, kd=qd)

        self.log('Cutoff energies:')
        self.log('    Wave functions:       %10.3f eV' %
                 (self.pd.ecut * Hartree))
        self.log('    Density:              %10.3f eV' %
                 (self.pd2.ecut * Hartree))

        # Calculate 1/|G+q|^2 with special treatment of |G+q|=0:
        G2_qG = self.pd2.G2_qG
        if self.q0 is None:
            if self.omega is None:
                self.iG2_qG = [1.0 / G2_G for G2_G in G2_qG]
            else:
                self.iG2_qG = [(1.0 / G2_G *
                                (1 - np.exp(-G2_G / (4 * self.omega**2))))
                               for G2_G in G2_qG]
        else:
            G2_qG[self.q0][0] = 117.0  # avoid division by zero
            if self.omega is None:
                self.iG2_qG = [1.0 / G2_G for G2_G in G2_qG]
                self.iG2_qG[self.q0][0] = self.gamma
            else:
                self.iG2_qG = [(1.0 / G2_G *
                                (1 - np.exp(-G2_G / (4 * self.omega**2))))
                               for G2_G in G2_qG]
                self.iG2_qG[self.q0][0] = 1 / (4 * self.omega**2)
            G2_qG[self.q0][0] = 0.0  # restore correct value

        # Compensation charges:
        self.ghat = PWLFC([setup.ghat_l for setup in wfs.setups], self.pd2)
        self.ghat.set_positions(self.spos_ac)

        if self.molecule:
            self.initialize_gaussian()
            self.log('Value of beta parameter: %.3f 1/Bohr^2' % self.beta)
            
        self.timer.stop('Initialization')
        
        # Ready ... set ... go:
        self.t0 = time()
        self.npairs = 0
        self.evv = 0.0
        self.evvacdf = 0.0
        for s in range(self.wfs.nspins):
            kpt1_q = [KPoint(self.wfs, noccmax).initialize(kpt)
                      for kpt in self.wfs.kpt_u if kpt.s == s]
            kpt2_q = kpt1_q[:]

            if len(kpt1_q) == 0:
                # No s-spins on this CPU:
                continue

            # Send and receive ranks:
            srank = self.wfs.kd.get_rank_and_index(
                s, (kpt1_q[0].k - 1) % kd.nibzkpts)[0]
            rrank = self.wfs.kd.get_rank_and_index(
                s, (kpt1_q[-1].k + 1) % kd.nibzkpts)[0]

            # Shift k-points kd.nibzkpts - 1 times:
            for i in range(kd.nibzkpts):
                if i < kd.nibzkpts - 1:
                    if kparallel:
                        kpt = kpt2_q[-1].next(self.wfs)
                        kpt.start_receiving(rrank)
                        kpt2_q[0].start_sending(srank)
                    else:
                        kpt = kpt2_q[0]

                self.timer.start('Calculate')
                for kpt1, kpt2 in zip(kpt1_q, kpt2_q):
                    # Loop over all k-points that k2 can be mapped to:
                    for K2, q, n1_n, n2 in self.indices(s, kpt1.k, kpt2.k):
                        self.apply(K2, q, kpt1, kpt2, n1_n, n2)
                self.timer.stop('Calculate')

                if i < kd.nibzkpts - 1:
                    self.timer.start('Wait')
                    if kparallel:
                        kpt.wait()
                        kpt2_q[0].wait()
                    self.timer.stop('Wait')
                    kpt2_q.pop(0)
                    kpt2_q.append(kpt)

        self.evv = self.world.sum(self.evv)
        self.evvacdf = self.world.sum(self.evvacdf)
        self.calculate_exx_paw_correction()
        
        if self.method == 'standard':
            self.exx = self.evv + self.devv + self.evc + self.ecc
        elif self.method == 'acdf':
            self.exx = self.evvacdf + self.devv + self.evc + self.ecc
        else:
            1 / 0

        self.log('Exact exchange energy:')
        for txt, e in [
            ('core-core', self.ecc),
            ('valence-core', self.evc),
            ('valence-valence (pseudo, acdf)', self.evvacdf),
            ('valence-valence (pseudo, standard)', self.evv),
            ('valence-valence (correction)', self.devv),
            ('total (%s)' % self.method, self.exx)]:
            self.log('    %-36s %14.6f eV' % (txt + ':', e * Hartree))

        self.log('Total time: %10.3f seconds' % (time() - self.t0))

        self.npairs = self.world.sum(self.npairs)
        assert self.npairs == self.npairs0
        
        self.timer.stop('EXX')
        self.timer.write(self.fd)

    def calculate_gamma(self, vol, alpha):
        if self.molecule:
            return 0.0

        N_c = self.kd.N_c
        offset_c = (N_c + 1) % 2 * 0.5 / N_c
        bzq_qc = monkhorst_pack(N_c) + offset_c
        qd = KPointDescriptor(bzq_qc)
        pd = PWDescriptor(self.wfs.pd.ecut, self.wfs.gd, kd=qd)
        gamma = (vol / (2 * pi)**2 * sqrt(pi / alpha) *
                 self.kd.nbzkpts)
        for G2_G in pd.G2_qG:
            if G2_G[0] < 1e-7:
                G2_G = G2_G[1:]
            gamma -= np.dot(np.exp(-alpha * G2_G), G2_G**-1)
        return gamma / self.qstride_c.prod()

    def indices(self, s, k1, k2):
        """Generator for (K2, q, n1, n2) indices for (k1, k2) pair.

        s: int
            Spin index.
        k1: int
            Index of k-point in the IBZ.
        k2: int
            Index of k-point in the IBZ.

        Returns (K, q, n1_n, n2), where K then index of the k-point in
        the BZ that k2 is mapped to, q is the index of the q-vector
        between K and k1, and n1_n is a list of bands that should be
        combined with band n2."""

        for K, k in enumerate(self.kd.bz2ibz_k):
            if k == k2:
                for K, q, n1_n, n2 in self._indices(s, k1, k2, K):
                    yield K, q, n1_n, n2
            
    def _indices(self, s, k1, k2, K2):
        k1_c = self.kd.ibzk_kc[k1]
        k2_c = self.kd.bzk_kc[K2]
        q_c = k2_c - k1_c
        q = abs(self.bzq_qc - q_c).sum(1).argmin()
        if abs(self.bzq_qc[q] - q_c).sum() > 1e-7:
            return

        if self.gamma_point == 0 and q == self.q0:
            return

        nocc1 = self.nocc_sk[s, k1]
        nocc2 = self.nocc_sk[s, k2]

        # Is k2 in the IBZ?
        is_ibz2 = (self.kd.ibz2bz_k[k2] == K2)

        for n2 in range(self.wfs.bd.nbands):
            # Find range of n1's (from n1a to n1b-1):
            if is_ibz2:
                # We get this combination twice, so let's only do half:
                if k1 >= k2:
                    n1a = n2
                else:
                    n1a = n2 + 1
            else:
                n1a = 0

            n1b = self.wfs.bd.nbands

            if self.bandstructure:
                if n2 >= nocc2:
                    n1b = min(n1b, nocc1)
            else:
                if n2 >= nocc2:
                    break
                n1b = min(n1b, nocc1)

            if self.bands is not None:
                assert self.bandstructure
                n1_n = []
                for n1 in range(n1a, n1b):
                    if (n1 in self.bands and n2 < nocc2 or
                        is_ibz2 and n2 in self.bands and n1 < nocc1):
                        n1_n.append(n1)
                n1_n = np.array(n1_n)
            else:
                n1_n = np.arange(n1a, n1b)

            if len(n1_n) == 0:
                continue

            yield K2, q, n1_n, n2

    def apply(self, K2, q, kpt1, kpt2, n1_n, n2):
        k20_c = self.kd.ibzk_kc[kpt2.k]
        k2_c = self.kd.bzk_kc[K2]

        if k2_c.any():
            self.timer.start('Initialize plane waves')
            eik2r_R = self.wfs.gd.plane_wave(k2_c)
            eik20r_R = self.wfs.gd.plane_wave(k20_c)
            self.timer.stop('Initialize plane waves')
        else:
            eik2r_R = 1.0
            eik20r_R = 1.0

        w1 = self.kd.weight_k[kpt1.k]
        w2 = self.kd.weight_k[kpt2.k]

        # Is k2 in the 1. BZ?
        is_ibz2 = (self.kd.ibz2bz_k[kpt2.k] == K2)

        e_n = self.calculate_interaction(n1_n, n2, kpt1, kpt2, q, K2,
                                         eik20r_R, eik2r_R,
                                         is_ibz2)

        e_n *= 1.0 / self.kd.nbzkpts / self.wfs.nspins * self.qstride_c.prod()
        
        if q == self.q0:
            e_n[n1_n == n2] *= 0.5

        f1_n = kpt1.f_n[n1_n]
        eps1_n = kpt1.eps_n[n1_n]
        f2 = kpt2.f_n[n2]
        eps2 = kpt2.eps_n[n2]

        s_n = np.sign(eps2 - eps1_n)

        evv = (f1_n * f2 * e_n).sum()
        evvacdf = 0.5 * (f1_n * (1 - s_n) * e_n +
                         f2 * (1 + s_n) * e_n).sum()
        self.evv += evv * w1
        self.evvacdf += evvacdf * w1
        if is_ibz2:
            self.evv += evv * w2
            self.evvacdf += evvacdf * w2

        if self.bandstructure:
            x = self.wfs.nspins
            self.exx_skn[kpt1.s, kpt1.k, n1_n] += x * f2 * e_n
            if is_ibz2:
                self.exx_skn[kpt2.s, kpt2.k, n2] += x * np.dot(f1_n, e_n)

    def calculate_interaction(self, n1_n, n2, kpt1, kpt2, q, k,
                              eik20r_R, eik2r_R, is_ibz2):
        """Calculate Coulomb interactions.

        For all n1 in the n1_n list, calculate interaction with n2."""

        # number of plane waves:
        ng1 = self.wfs.ng_k[kpt1.k]
        ng2 = self.wfs.ng_k[kpt2.k]

        # Transform to real space and apply symmetry operation:
        self.timer.start('IFFT1')
        if is_ibz2:
            u2_R = self.pd.ifft(kpt2.psit_nG[n2, :ng2], kpt2.k)
        else:
            psit2_R = self.pd.ifft(kpt2.psit_nG[n2, :ng2], kpt2.k) * eik20r_R
            self.timer.start('Symmetry transform')
            u2_R = self.kd.transform_wave_function(psit2_R, k) / eik2r_R
            self.timer.stop()
        self.timer.stop()

        # Calculate pair densities:
        nt_nG = self.pd2.zeros(len(n1_n), q=q)
        for n1, nt_G in zip(n1_n, nt_nG):
            self.timer.start('IFFT2')
            u1_R = self.pd.ifft(kpt1.psit_nG[n1, :ng1], kpt1.k)
            self.timer.stop()
            nt_R = u1_R.conj() * u2_R
            self.timer.start('FFT')
            nt_G[:] = self.pd2.fft(nt_R, q)
            self.timer.stop()
        
        s = self.kd.sym_k[k]
        time_reversal = self.kd.time_reversal_k[k]
        k2_c = self.kd.ibzk_kc[kpt2.k]

        self.timer.start('Compensation charges')
        Q_anL = {}  # coefficients for shape functions
        for a, P1_ni in kpt1.P_ani.items():
            P1_ni = P1_ni[n1_n]

            if is_ibz2:
                P2_i = kpt2.P_ani[a][n2]
            else:
                b = self.kd.symmetry.a_sa[s, a]
                S_c = (np.dot(self.spos_ac[a], self.kd.symmetry.op_scc[s]) -
                       self.spos_ac[b])
                assert abs(S_c.round() - S_c).max() < 1e-5
                if self.ghat.dtype == complex:
                    x = np.exp(2j * pi * np.dot(k2_c, S_c))
                else:
                    x = 1.0
                P2_i = np.dot(self.wfs.setups[a].R_sii[s],
                              kpt2.P_ani[b][n2]) * x
                if time_reversal:
                    P2_i = P2_i.conj()

            D_np = []
            for P1_i in P1_ni:
                D_ii = np.outer(P1_i.conj(), P2_i)
                D_np.append(pack(D_ii))
            Q_anL[a] = np.dot(D_np, self.wfs.setups[a].Delta_pL)
            
        self.timer.start('Expand')
        if q != self.qlatest:
            self.f_IG = self.ghat.expand(q)
            self.qlatest = q
        self.timer.stop('Expand')

        # Add compensation charges:
        self.ghat.add(nt_nG, Q_anL, q, self.f_IG)
        self.timer.stop('Compensation charges')

        if self.molecule and n2 in n1_n:
            nn = (n1_n == n2).nonzero()[0][0]
            nt_nG[nn] -= self.ngauss_G
        else:
            nn = None
            
        iG2_G = self.iG2_qG[q]
        
        # Calculate energies:
        e_n = np.empty(len(n1_n))
        for n, nt_G in enumerate(nt_nG):
            e_n[n] = -4 * pi * np.real(self.pd2.integrate(nt_G, nt_G * iG2_G))
            self.npairs += 1
        
        if nn is not None:
            e_n[nn] -= 2 * (self.pd2.integrate(nt_nG[nn], self.vgauss_G) +
                            (self.beta / 2 / pi)**0.5)

        if self.write_timing_information:
            t = (time() - self.t0) / len(n1_n)
            self.log('Time for first pair-density: %10.3f seconds' % t)
            self.log('Estimated total time:        %10.3f seconds' %
                     (t * self.npairs0 / self.world.size))
            self.write_timing_information = False

        return e_n

    def calculate_exx_paw_correction(self):
        self.timer.start('PAW correction')
        self.devv = 0.0
        self.evc = 0.0
        self.ecc = 0.0
                         
        deg = 2 // self.wfs.nspins  # spin degeneracy
        for a, D_sp in self.dens.D_asp.items():
            setup = self.wfs.setups[a]
            for D_p in D_sp:
                D_ii = unpack2(D_p)
                ni = len(D_ii)

                for i1 in range(ni):
                    for i2 in range(ni):
                        A = 0.0
                        for i3 in range(ni):
                            p13 = packed_index(i1, i3, ni)
                            for i4 in range(ni):
                                p24 = packed_index(i2, i4, ni)
                                A += setup.M_pp[p13, p24] * D_ii[i3, i4]
                        self.devv -= D_ii[i1, i2] * A / deg

                self.evc -= np.dot(D_p, setup.X_p)
            self.ecc += setup.ExxC

        if not self.bandstructure:
            self.timer.stop('PAW correction')
            return

        Q = self.world.size // self.wfs.kd.comm.size
        self.exx_skn *= Q
        for kpt in self.wfs.kpt_u:
            for a, D_sp in self.dens.D_asp.items():
                setup = self.wfs.setups[a]
                for D_p in D_sp:
                    D_ii = unpack2(D_p)
                    ni = len(D_ii)
                    P_ni = kpt.P_ani[a]
                    for i1 in range(ni):
                        for i2 in range(ni):
                            A = 0.0
                            for i3 in range(ni):
                                p13 = packed_index(i1, i3, ni)
                                for i4 in range(ni):
                                    p24 = packed_index(i2, i4, ni)
                                    A += setup.M_pp[p13, p24] * D_ii[i3, i4]
                            self.exx_skn[kpt.s, kpt.k] -= \
                                (A * P_ni[:, i1].conj() * P_ni[:, i2]).real
                            p12 = packed_index(i1, i2, ni)
                            self.exx_skn[kpt.s, kpt.k] -= \
                                (P_ni[:, i1].conj() * setup.X_p[p12] *
                                 P_ni[:, i2]).real / self.wfs.nspins

        self.world.sum(self.exx_skn)
        self.exx_skn *= self.hybrid / Q
        self.timer.stop('PAW correction')
    
    def initialize_gaussian(self):
        """Calculate gaussian compensation charge and its potential.

        Used to decouple electrostatic interactions between
        periodically repeated images for molecular calculations.

        Charge containing one electron::

            (beta/pi)^(3/2)*exp(-beta*r^2),

        its Fourier transform::

            exp(-G^2/(4*beta)),

        and its potential::

            erf(beta^0.5*r)/r.
        """

        gd = self.wfs.gd

        # Set exponent of exp-function to -19 on the boundary:
        self.beta = 4 * 19 * (gd.icell_cv**2).sum(1).max()

        # Calculate gaussian:
        G_Gv = self.pd2.get_reciprocal_vectors()
        G2_G = self.pd2.G2_qG[0]
        C_v = gd.cell_cv.sum(0) / 2  # center of cell
        self.ngauss_G = np.exp(-1.0 / (4 * self.beta) * G2_G +
                                1j * np.dot(G_Gv, C_v)) / gd.dv

        # Calculate potential from gaussian:
        R_Rv = gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        r_R = ((R_Rv - C_v)**2).sum(3)**0.5
        if (gd.N_c % 2 == 0).all():
            r_R[tuple(gd.N_c // 2)] = 1.0  # avoid dividing by zero
        v_R = erf(self.beta**0.5 * r_R) / r_R
        if (gd.N_c % 2 == 0).all():
            v_R[tuple(gd.N_c // 2)] = (4 * self.beta / pi)**0.5
        self.vgauss_G = self.pd2.fft(v_R)

        # Compare self-interaction to analytic result:
        assert abs(0.5 * self.pd2.integrate(self.ngauss_G, self.vgauss_G) -
                   (self.beta / 2 / pi)**0.5) < 1e-6


class KPoint:
    def __init__(self, wfs, nbands):
        """Helper class for parallelizing over k-points.

        Placeholder for wave functions, occupation numbers, eigenvalues,
        projections, spin index and global k-point index."""
        
        self.kd = wfs.kd
        self.ng_k = wfs.ng_k

        # Array large enough to hold wave functions from all k-points:
        self.psit_nG = wfs.pd.empty(nbands)

        self.requests = []

    def initialize(self, kpt):
        ng = self.ng_k[kpt.k]
        nbands = len(self.psit_nG)
        self.psit_nG[:, :ng] = kpt.psit_nG[:nbands]
        self.f_n = kpt.f_n / kpt.weight  # will be in the range [0,1]
        self.eps_n = kpt.eps_n
        self.P_ani = kpt.P_ani
        self.k = kpt.k
        self.s = kpt.s

        return self

    def next(self, wfs):
        """Create empty object.

        Data will be received from another process."""
        
        nbands = len(self.psit_nG)
        kpt = KPoint(wfs, nbands)

        # Allocate arrays for receiving:
        kpt.f_n = wfs.bd.empty()
        kpt.eps_n = wfs.bd.empty()

        # Total number of projector functions:
        I = sum([P_ni.shape[1] for P_ni in self.P_ani.values()])
        
        kpt.P_nI = np.empty((wfs.bd.nbands, I), wfs.dtype)

        kpt.P_ani = {}
        I1 = 0
        assert self.P_ani.keys() == range(len(self.P_ani))  # ???
        for a, P_ni in self.P_ani.items():
            I2 = I1 + P_ni.shape[1]
            kpt.P_ani[a] = kpt.P_nI[:, I1:I2]
            I1 = I2

        kpt.k = (self.k + 1) % self.kd.nibzkpts
        kpt.s = self.s
        
        return kpt
        
    def start_sending(self, rank):
        assert self.P_ani.keys() == range(len(self.P_ani))  # ???
        P_nI = np.hstack([P_ni for P_ni in self.P_ani.values()])
        P_nI = np.ascontiguousarray(P_nI)
        self.requests += [
            self.kd.comm.send(self.psit_nG, rank, block=False, tag=1),
            self.kd.comm.send(self.f_n, rank, block=False, tag=2),
            self.kd.comm.send(self.eps_n, rank, block=False, tag=3),
            self.kd.comm.send(P_nI, rank, block=False, tag=4)]
        
    def start_receiving(self, rank):
        self.requests += [
            self.kd.comm.receive(self.psit_nG, rank, block=False, tag=1),
            self.kd.comm.receive(self.f_n, rank, block=False, tag=2),
            self.kd.comm.receive(self.eps_n, rank, block=False, tag=3),
            self.kd.comm.receive(self.P_nI, rank, block=False, tag=4)]
    
    def wait(self):
        self.kd.comm.waitall(self.requests)
        self.requests = []
