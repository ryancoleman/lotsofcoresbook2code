# Copyright (C) 2010  CAMd
# Please see the accompanying LICENSE file for further information.

"""This module provides all the classes and functions associated with the
evaluation of exact exchange with k-point sampling."""

from math import pi, sqrt
from time import time

import numpy as np
from ase.utils import prnt

from gpaw.utilities import pack, unpack2, packed_index, logfile
from gpaw.lfc import LFC
from gpaw.wavefunctions.pw import PWDescriptor, PWWaveFunctions
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.xc.hybrid import HybridXCBase

class KPoint:
    def __init__(self, kd, kpt=None):
        """Helper class for parallelizing over k-points.

        Placeholder for wave functions, occupation numbers,
        projections, and global k-point index."""
        
        self.kd = kd
        
        if kpt is not None:
            self.psit_nG = kpt.psit_nG
            self.f_n = kpt.f_n / kpt.weight / kd.nbzkpts * 2 / kd.nspins
            self.eps_n = kpt.eps_n
            self.weight = np.array([kpt.weight])

            self.P_ani = kpt.P_ani
            self.k = kpt.k
            self.s = kpt.s
 
        self.requests = []
        
    def next(self):
        """Create empty object.

        Data will be received from other processor."""
        
        kpt = KPoint(self.kd)

        # intialize array for receiving:
        kpt.psit_nG = np.empty_like(self.psit_nG)
        kpt.f_n = np.empty_like(self.f_n)
        kpt.eps_n = np.empty_like(self.eps_n)
        kpt.weight = np.empty_like(self.weight)

        # Total number of projector functions:
        I = sum([P_ni.shape[1] for P_ni in self.P_ani.values()])
        
        kpt.P_nI = np.empty((len(kpt.f_n), I), complex)

        kpt.P_ani = {}
        I1 = 0
        for a, P_ni in self.P_ani.items():
            I2 = I1 + P_ni.shape[1]
            kpt.P_ani[a] = kpt.P_nI[:, I1:I2]
            I1 = I2

        kpt.k = (self.k + 1) % self.kd.nibzkpts
        kpt.s = self.s
        
        return kpt
        
    def start_sending(self, rank):
        P_nI = np.hstack([P_ni for P_ni in self.P_ani.values()])
        P_nI = np.ascontiguousarray(P_nI)
        self.requests += [
            self.kd.comm.send(self.psit_nG, rank, block=False, tag=1),
            self.kd.comm.send(self.f_n, rank, block=False, tag=2),
            self.kd.comm.send(self.eps_n, rank, block=False, tag=3),
            self.kd.comm.send(self.weight, rank, block=False, tag=4),
            self.kd.comm.send(P_nI, rank, block=False, tag=5)]
        
    def start_receiving(self, rank):
        self.requests += [
            self.kd.comm.receive(self.psit_nG, rank, block=False, tag=1),
            self.kd.comm.receive(self.f_n, rank, block=False, tag=2),
            self.kd.comm.receive(self.eps_n, rank, block=False, tag=3),
            self.kd.comm.receive(self.weight, rank, block=False, tag=4),
            self.kd.comm.receive(self.P_nI, rank, block=False, tag=5)]

        
    def wait(self):
        self.kd.comm.waitall(self.requests)
        self.requests = []
        

class HybridXC(HybridXCBase):
    orbital_dependent = True
    def __init__(self, name, hybrid=None, xc=None, gygi=False,
                 alpha=None, skip_gamma=False, ecut=None, 
                 etotflag = False, acdf=False, coredensity=True,
                 logfilename='-', bands=None, core_valence=True):
        """Mix standard functionals with exact exchange.

        bands: list or None
            List of bands to calculate energy for.  Default is None
            meaning do all bands.
        """

        self.alpha = alpha
        self.skip_gamma = skip_gamma
        self.gygi = gygi
        
        self.exx = 0.0
        self.etotflag = etotflag
        self.ecut = ecut
        self.fd = logfilename
        self.write_timing_information = True
        self.bands = bands
        self.acdf = acdf  # adiabatic-connection dissipation fluctuation for RPA correlation energy
        self.coredensity = coredensity
        self.core_valence = core_valence
        if self.acdf:
            self.exxacdf = 0.0
            self.etotflag = True
            print('etotflag is True')

        HybridXCBase.__init__(self, name, hybrid, xc)

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
        addcoredensity = self.coredensity # XXX overwrites input

        return self.xc.calculate_paw_correction(setup, D_sp, dEdD_sp,
                                 addcoredensity, a)
    
    def initialize(self, density, hamiltonian, wfs, occupations):
        self.xc.initialize(density, hamiltonian, wfs, occupations)
        self.nspins = wfs.nspins
        self.setups = wfs.setups
        self.density = density
        self.kpt_u = wfs.kpt_u
        
        self.gd = density.gd
        self.kd = wfs.kd
        self.bd = wfs.bd
        if self.bd.comm.size > 1:
            raise ValueError('Band parallelization not supported by hybridk')
        self.wfs = wfs

        self.world = wfs.world

        self.fd = logfile(self.fd, self.world.rank)

        N = self.gd.N_c.prod()
        vol = self.gd.dv * N
        
        if self.alpha is None:
            # XXX ?
            self.alpha = 6 * vol**(2 / 3.0) / pi**2
            
        if self.ecut is None:
            self.ecut = 0.5 * pi**2 / (self.gd.h_cv**2).sum(1).max() * 0.9999
            
        self.bzq_qc = self.kd.get_bz_q_points()
        qd = KPointDescriptor(self.bzq_qc)
        q0 = self.kd.where_is_q(np.zeros(3), self.bzq_qc)
        
        self.pwd = PWDescriptor(self.ecut, self.gd, complex, kd=qd)

        G2_qG = self.pwd.G2_qG
        G2_qG[q0][0] = 117.0
        self.iG2_qG = [1.0 / G2_G for G2_G in G2_qG]
        G2_qG[q0][0] = 0.0
        self.iG2_qG[q0][0] = 0.0

        self.gamma = (vol / (2 * pi)**2 * sqrt(pi / self.alpha) *
                      self.kd.nbzkpts)

        for q in range(self.kd.nbzkpts):
            self.gamma -= np.dot(np.exp(-self.alpha * G2_qG[q]),
                                 self.iG2_qG[q])

        self.iG2_qG[q0][0] = self.gamma
        
        self.ghat = LFC(self.gd,
                        [setup.ghat_l for setup in density.setups],
                        qd, dtype=complex)
        
        self.log('Value of alpha parameter:', self.alpha)
        self.log('Value of gamma parameter:', self.gamma)
        self.log('Cutoff energy:', self.ecut, 'Hartree')
        self.log('%d x %d x %d k-points' % tuple(self.kd.N_c))

    def set_positions(self, spos_ac):
        self.ghat.set_positions(spos_ac)
        self.spos_ac = spos_ac

    def calculate(self, gd, n_sg, v_sg=None, e_g=None):
        # Normal XC contribution:
        exc = self.xc.calculate(gd, n_sg, v_sg, e_g)

         # Add EXX contribution:
        return exc + self.exx

    def calculate_exx(self):
        """Non-selfconsistent calculation."""

        kd = self.kd
        K = kd.nibzkpts
        W = self.world.size // self.nspins
        parallel = (W > 1)

        self.log("%d CPU's used for %d IBZ k-points" % (W, K))
        self.log('Spins:', self.nspins)

        if self.etotflag and not self.gygi:
            self.nbandstmp = 0
            for s in range(self.nspins):
                kpt1_k = [KPoint(kd, kpt)
                          for kpt in self.kpt_u if kpt.s == s]
                for kpt1 in kpt1_k:
                    for n1 in range(self.bd.nbands):
                        f_n = kpt1.f_n[n1]
                        if np.abs(f_n) < 1e-10:
                            self.nbandstmp = max(self.nbandstmp, n1)
                            break
                    else:
                        self.nbandstmp = self.bd.nbands

            tmp = np.zeros(kd.comm.size, dtype=int)
            kd.comm.all_gather(np.array([self.nbandstmp]), tmp)
            self.nbands = tmp.max()
        else:
            self.nbands = self.bd.nbands
        
        B = self.nbands
        self.log('Number of bands calculated:', B)
        self.log('Number of valence electrons:', self.setups.nvalence)

        E = B - self.setups.nvalence / 2.0  # empty bands
        self.npairs = (K * kd.nbzkpts - 0.5 * K**2) * (B**2 - E**2)
        self.log('Approximate number of pairs:', self.npairs)

        if not self.etotflag:
            self.exx_skn = np.zeros((self.nspins, K, B))
            self.debug_skn = np.zeros((self.nspins, K, B))

        for s in range(self.nspins):
            kpt1_q = [KPoint(kd, kpt)
                      for kpt in self.kpt_u if kpt.s == s]
            kpt2_q = kpt1_q[:]

            if len(kpt1_q) == 0:
                # No s-spins on this CPU:
                continue

            # Send rank:
            srank = kd.get_rank_and_index(s, (kpt1_q[0].k - 1) % K)[0]
            # Receive rank:
            rrank = kd.get_rank_and_index(s, (kpt1_q[-1].k + 1) % K)[0]

            # Shift k-points K - 1 times:
            for i in range(K):
                if i < K - 1:
                    if parallel:
                        kpt = kpt2_q[-1].next()
                        kpt.start_receiving(rrank)
                        kpt2_q[0].start_sending(srank)
                    else:
                        kpt = kpt2_q[0]

                for kpt1, kpt2 in zip(kpt1_q, kpt2_q):
                    for k, ik in enumerate(kd.bz2ibz_k):
                        if ik == kpt2.k:
                            self.apply(kpt1, kpt2, k)

                if i < K - 1:
                    if parallel:
                        kpt.wait()
                        kpt2_q[0].wait()
                    kpt2_q.pop(0)
                    kpt2_q.append(kpt)

        if self.etotflag:
            if self.acdf:
                self.exxacdf = self.world.sum(self.exxacdf[0])
                self.exx = self.exxacdf
            else:
                self.exx = self.world.sum(self.exx)
            self.exx += self.calculate_exx_paw_correction()
                
        else:
            for kpt in self.kpt_u:
                for a, D_sp in self.density.D_asp.items():
                    setup = self.setups[a]
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
                                    (self.hybrid * A *
                                     P_ni[:, i1].conj() * P_ni[:, i2]).real
                                
                                p12 = packed_index(i1, i2, ni)
                                if self.core_valence:
                                    if setup.X_p is not None:
                                        self.exx_skn[kpt.s, kpt.k] -= self.hybrid * \
                                                                      (P_ni[:, i1].conj() * setup.X_p[p12] *
                                                                       P_ni[:, i2]).real / self.nspins

            
            self.world.sum(self.exx_skn)
            self.exx = 0.0
            for kpt in self.kpt_u:
                self.exx += 0.5 * np.dot(kpt.f_n, self.exx_skn[kpt.s, kpt.k])
            self.exx = self.world.sum(self.exx)

            for a, D_sp in self.density.D_asp.items():
                setup = self.setups[a]
                if self.coredensity:
                    self.exx += self.hybrid * setup.ExxC
                if self.core_valence:
                    self.exx -= self.hybrid * 0.5 * np.dot(D_sp.sum(0), setup.X_p)

            self.world.sum(self.debug_skn)
            assert (self.debug_skn == self.kd.nbzkpts * B).all()
    
    def apply(self, kpt1, kpt2, k):
        k1_c = self.kd.ibzk_kc[kpt1.k]
        k20_c = self.kd.ibzk_kc[kpt2.k]
        k2_c = self.kd.bzk_kc[k]
        q_c = k2_c - k1_c
        N_c = self.gd.N_c

        q = self.kd.where_is_q(q_c, self.bzq_qc)
        
        q_c = self.bzq_qc[q]
        eik1r_R = np.exp(2j * pi * np.dot(np.indices(N_c).T, k1_c / N_c).T)
        eik2r_R = np.exp(2j * pi * np.dot(np.indices(N_c).T, k20_c / N_c).T)
        eiqr_R = np.exp(2j * pi * np.dot(np.indices(N_c).T, q_c / N_c).T)

        same = abs(k1_c - k2_c).max() < 1e-9

        iG2_G = self.iG2_qG[q]
            
        N = N_c.prod()
        vol = self.gd.dv * N
        nspins = self.nspins

        fcut = 1e-10
        is_ibz2 = abs(k2_c - self.kd.ibzk_kc[kpt2.k]).max() < 1e-9
        
        for n1 in range(self.nbands):
            f1 = kpt1.f_n[n1]
            e1 = kpt1.eps_n[n1]
            for n2 in range(self.nbands):
                if same:
                    assert is_ibz2
                    if n2 > n1:
                        continue
                elif is_ibz2:
                    if kpt1.k > kpt2.k:
                        if n2 > n1:
                            continue
                    else:
                        if n2 >= n1:
                            continue
                        
                f2 = kpt2.f_n[n2]
                e2 = kpt2.eps_n[n2]

                x = 1.0
                if same and n1 == n2:
                    x = 0.5

                if not self.etotflag:
                    self.debug_skn[kpt1.s, kpt1.k, n1] += x
                    if is_ibz2:
                        self.debug_skn[kpt2.s, kpt2.k, n2] += x

                if self.etotflag and not self.gygi:
                    if abs(f1) < fcut or abs(f2) < fcut:
                        continue
                else:
                    if abs(f1) < fcut and abs(f2) < fcut:
                        continue

                if self.bands is not None:
                    if not (n1 in self.bands or is_ibz2 and n2 in self.bands):
                        continue

                if self.skip_gamma and same:
                    continue
                
                t0 = time()
                nt_R = self.calculate_pair_density(n1, n2, kpt1, kpt2, q, k,
                                                   eik1r_R, eik2r_R, eiqr_R,
                                                   is_ibz2)
                nt_G = self.pwd.fft(nt_R, q) / N
                vt_G = nt_G.copy()
                vt_G *= -pi * vol * iG2_G
                e = np.vdot(nt_G, vt_G).real * nspins * self.hybrid * x
                
                if self.etotflag:
                    if self.acdf:
                        if self.gygi and same:
                            self.exxacdf += f2 * e * kpt1.weight                            
                        else:
                            self.exxacdf += 0.5 * (f1 * (1-np.sign(e2-e1)) * e + 
                                                   f2 * (1-np.sign(e1-e2)) * e ) * kpt1.weight
                    else:
                        self.exx += f2 * e * kpt1.weight[0] * f1 * self.kd.nbzkpts * nspins / 2
                else:
                    self.exx_skn[kpt1.s, kpt1.k, n1] += 2 * f2 * e

                if is_ibz2:
                    if self.etotflag:
                        if self.acdf:
                            if self.gygi and same:
                                self.exxacdf += f1 * e * kpt2.weight
                            else:
                                self.exxacdf += 0.5 * (f1 * (1-np.sign(e2-e1)) * e +
                                                       f2 * (1-np.sign(e1-e2)) * e ) * kpt2.weight
                        else:
                            self.exx += f1 * e * kpt2.weight[0] * f2 * self.kd.nbzkpts * nspins / 2
                    else:
                        self.exx_skn[kpt2.s, kpt2.k, n2] += 2 * f1 * e
                    
                if self.write_timing_information:
                    t = time() - t0
                    self.log('Time for first pair-density:', t, 'seconds')
                    self.log('Estimated total time',
                             t * self.npairs / self.world.size, 'seconds')
                    self.write_timing_information = False

    def calculate_exx_paw_correction(self):
        exx = 0
        deg = 2 // self.nspins  # spin degeneracy
        for a, D_sp in self.density.D_asp.items():
            setup = self.setups[a]
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
                        exx -= self.hybrid / deg * D_ii[i1, i2] * A

                if self.core_valence:
                    if setup.X_p is not None:
                        exx -= self.hybrid * np.dot(D_p, setup.X_p)
            if self.coredensity:
                exx += self.hybrid * setup.ExxC
        return exx

    def calculate_pair_density(self, n1, n2, kpt1, kpt2, q, k,
                               eik1r_R, eik2r_R, eiqr_R, ibz2):
        if isinstance(self.wfs, PWWaveFunctions):
            psit1_R = self.wfs.pd.ifft(kpt1.psit_nG[n1]) * eik1r_R
            psit2_R = self.wfs.pd.ifft(kpt2.psit_nG[n2]) * eik2r_R
        else:
            psit1_R = kpt1.psit_nG[n1]
            psit2_R = kpt2.psit_nG[n2]

        if ibz2:
            psit2_R = psit2_R
        else:
            psit2_R = np.asarray(self.kd.transform_wave_function(psit2_R, k),
                                 complex)
        nt_R = psit1_R.conj() * psit2_R

        s = self.kd.sym_k[k]
        time_reversal = self.kd.time_reversal_k[k]
        k2_c = self.kd.ibzk_kc[kpt2.k]

        Q_aL = {}
        for a, P1_ni in kpt1.P_ani.items():
            P1_i = P1_ni[n1]

            b = self.kd.symmetry.a_sa[s, a]
            S_c = (np.dot(self.spos_ac[a], self.kd.symmetry.op_scc[s]) -
                   self.spos_ac[b])
            assert abs(S_c.round() - S_c).max() < 1e-13
            x = np.exp(2j * pi * np.dot(k2_c, S_c))
            P2_i = np.dot(self.setups[a].R_sii[s], kpt2.P_ani[b][n2]) * x
            if time_reversal:
                P2_i = P2_i.conj()

            if ibz2:
                P2_i = kpt2.P_ani[a][n2]

            D_ii = np.outer(P1_i.conj(), P2_i)
            D_p = pack(D_ii)
            Q_aL[a] = np.dot(D_p, self.setups[a].Delta_pL)

        self.ghat.add(nt_R, Q_aL, q)
        return nt_R / eiqr_R
