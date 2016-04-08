# Copyright (C) 2010  CAMd
# Please see the accompanying LICENSE file for further information.

"""This module provides all the classes and functions associated with the
evaluation of exact exchange with k-point sampling."""
from __future__ import print_function
from math import pi, sqrt
import sys
from time import ctime
import numpy as np

from ase import Atoms
from ase.units import Ha
from ase.utils import devnull

from gpaw.xc import XC
from gpaw.xc.kernel import XCNull
from gpaw.xc.functional import XCFunctional
from gpaw.utilities import pack, unpack2, packed_index
from gpaw.lfc import LFC
from gpaw.wavefunctions.pw import PWDescriptor
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.mpi import world, rank

class KPoint:
    def __init__(self, kd, kpt=None):
        """Helper class for parallelizing over k-points.

        Placeholder for wave functions, occupation numbers,
        projections, and global k-point index."""
        
        self.kd = kd
        
        if kpt is not None:
            
            self.psit_nG = kpt.psit_nG
            self.f_n = kpt.f_n / kpt.weight / kd.nbzkpts * 2 / kd.nspins
            self.weight = 1. / kd.nbzkpts * 2 / kd.nspins
            self.eps_n = kpt.eps_n
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

        # Total number of projector functions:
        I = sum([P_ni.shape[1] for P_ni in self.P_ani.values()])
        
        kpt.P_In = np.empty((I, len(kpt.f_n)), complex)

        kpt.P_ani = {}
        I1 = 0
        for a, P_ni in self.P_ani.items():
            I2 = I1 + P_ni.shape[1]
            kpt.P_ani[a] = kpt.P_In[I1:I2].T
            I1 = I2

        kpt.k = (self.k + 1) % self.kd.nibzkpts
        kpt.s = self.s
        
        return kpt
        
    def start_sending(self, rank):
        P_In = np.concatenate([P_ni.T for P_ni in self.P_ani.values()])
        self.requests += [
            self.kd.comm.send(self.psit_nG, rank, block=False, tag=1),
            self.kd.comm.send(self.f_n, rank, block=False, tag=2),
            self.kd.comm.send(P_In, rank, block=False, tag=3)]
        
    def start_receiving(self, rank):
        self.requests += [
            self.kd.comm.receive(self.psit_nG, rank, block=False, tag=1),
            self.kd.comm.receive(self.f_n, rank, block=False, tag=2),
            self.kd.comm.receive(self.P_In, rank, block=False, tag=3)]
        
    def wait(self):
        self.kd.comm.waitall(self.requests)
        self.requests = []
        

class HybridXC(XCFunctional):
    orbital_dependent = True
    def __init__(self, name, hybrid=None, xc=None, finegrid=False,
                 alpha=None, skip_gamma=False, gygi=False, acdf=True,
                 qsym=True, txt=None, ecut=None):
        """Mix standard functionals with exact exchange.

        name: str
            Name of hybrid functional.
        hybrid: float
            Fraction of exact exchange.
        xc: str or XCFunctional object
            Standard DFT functional with scaled down exchange.
        finegrid: boolean
            Use fine grid for energy functional evaluations?
        """

        if name == 'EXX':
            assert hybrid is None and xc is None
            hybrid = 1.0
            xc = XC(XCNull())
        elif name == 'PBE0':
            assert hybrid is None and xc is None
            hybrid = 0.25
            xc = XC('HYB_GGA_XC_PBEH')
        elif name == 'B3LYP':
            assert hybrid is None and xc is None
            hybrid = 0.2
            xc = XC('HYB_GGA_XC_B3LYP')
            
        if isinstance(xc, str):
            xc = XC(xc)

        self.hybrid = hybrid
        self.xc = xc
        self.type = xc.type
        self.alpha = alpha
        self.qsym = qsym
        self.skip_gamma = skip_gamma
        self.gygi = gygi
        self.acdf = acdf
        self.exx = None
        self.ecut = ecut
        if txt is None:
            if rank == 0:
                #self.txt = devnull
                self.txt = sys.stdout
            else:
                sys.stdout = devnull
                self.txt = devnull
        else:
            assert type(txt) is str
            from ase.parallel import paropen
            self.txt = paropen(txt, 'w')

        XCFunctional.__init__(self, name)
        
    def get_setup_name(self):
        return 'PBE'

    def calculate_radial(self, rgd, n_sLg, Y_L, v_sg,
                         dndr_sLg=None, rnablaY_Lv=None,
                         tau_sg=None, dedtau_sg=None):
        return self.xc.calculate_radial(rgd, n_sLg, Y_L, v_sg,
                                        dndr_sLg, rnablaY_Lv)

    def calculate_paw_correction(self, setup, D_sp, dEdD_sp=None,
                                 addcoredensity=True, a=None):
        return self.xc.calculate_paw_correction(setup, D_sp, dEdD_sp,
                                 addcoredensity, a)
    
    def initialize(self, density, hamiltonian, wfs, occupations):
        self.xc.initialize(density, hamiltonian, wfs, occupations)
        self.nspins = wfs.nspins
        self.setups = wfs.setups
        self.density = density
        self.kpt_u = wfs.kpt_u
        self.wfs = wfs
        
        self.gd = density.gd
        self.kd = wfs.kd
        self.bd = wfs.bd

        N_c = self.gd.N_c
        N = self.gd.N_c.prod()
        vol = self.gd.dv * N
        
        if self.alpha is None:
            # XXX ?
            self.alpha = 6 * vol**(2 / 3.0) / pi**2
            
        self.gamma = (vol / (2 * pi)**2 * sqrt(pi / self.alpha) *
                      self.kd.nbzkpts)

        if self.ecut is None:
            self.ecut = 0.5 * pi**2 / (self.gd.h_cv**2).sum(1).max() * 0.9999

        assert self.kd.N_c is not None
        n = self.kd.N_c * 2 - 1
        bzk_kc = np.indices(n).transpose((1, 2, 3, 0))
        bzk_kc.shape = (-1, 3)
        bzk_kc -= self.kd.N_c - 1
        self.bzk_kc = bzk_kc.astype(float) / self.kd.N_c

        self.bzq_qc = self.kd.get_bz_q_points()
        if self.qsym:
            op_scc = self.kd.symmetry.op_scc
            self.ibzq_qc = self.kd.get_ibz_q_points(self.bzq_qc,
                                                    op_scc)[0]
            self.q_weights = self.kd.q_weights * len(self.bzq_qc)
        else:
            self.ibzq_qc = self.bzq_qc
            self.q_weights = np.ones(len(self.bzq_qc))

        self.pwd = PWDescriptor(self.ecut, self.gd, complex)
        self.G2_qG = self.pwd.g2(self.bzk_kc)
        
        n = 0
        for k_c, Gpk2_G in zip(self.bzk_kc[:], self.G2_qG):
            if (k_c > -0.5).all() and (k_c <= 0.5).all(): #XXX???
                if k_c.any():
                    self.gamma -= np.dot(np.exp(-self.alpha * Gpk2_G),
                                         Gpk2_G**-1)
                else:
                    self.gamma -= np.dot(np.exp(-self.alpha * Gpk2_G[1:]),
                                         Gpk2_G[1:]**-1)
                n += 1
                
        assert n == self.kd.N_c.prod()

        self.pwd = PWDescriptor(self.ecut, self.gd, complex)
        self.G2_qG = self.pwd.g2(self.ibzq_qc)
        
        self.ghat = LFC(self.gd,
                        [setup.ghat_l for setup in density.setups],
                        KPointDescriptor(self.bzq_qc), dtype=complex)
        
        #self.interpolator = density.interpolator
        self.print_initialization(hamiltonian.xc.name)

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
        K = len(kd.bzk_kc)
        W = world.size // self.nspins
        parallel = (W > 1)
        
        self.exx = 0.0
        self.exx_kq = np.zeros((K, len(self.ibzq_qc)), float)
                
        for s in range(self.nspins):
            ibz_kpts = [KPoint(kd, kpt)
                        for kpt in self.kpt_u if kpt.s == s]
            for ik, kpt in enumerate(kd.bzk_kc):
                print('K %s %s ...' % (ik, kpt), file=self.txt)
                for iq, q in enumerate(self.ibzq_qc):
                    kpq = kd.find_k_plus_q(q, kpts_k=[ik])
                    self.apply(ibz_kpts[kd.bz2ibz_k[ik]],
                               ibz_kpts[kd.bz2ibz_k[kpq[0]]],
                               ik, kpq[0], iq)
                    
        self.exx = world.sum(self.exx)
        self.exx += self.calculate_exx_paw_correction()

        exx_q = np.sum(self.exx_kq, 0)

        print(file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print(file=self.txt)
        print('Contributions: q         w        E_q (eV)', file=self.txt) 
        for q in range(len(exx_q)):
            print('[%1.3f %1.3f %1.3f]    %1.3f   %s' % \
                  (self.ibzq_qc[q][0], self.ibzq_qc[q][1], self.ibzq_qc[q][2],
                   self.q_weights[q]/len(self.bzq_qc),
                   exx_q[q]/self.q_weights[q]*len(self.bzq_qc)*Ha), file=self.txt)
        print('E_EXX = %s eV' % (self.exx*Ha), file=self.txt)
        print(file=self.txt)
        print('Calculation completed at:  ', ctime(), file=self.txt)
        print(file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print(file=self.txt)
         
    def apply(self, kpt1, kpt2, ik1, ik2, iq):
        k1_c = self.kd.bzk_kc[ik1]
        k2_c = self.kd.bzk_kc[ik2]
        q = self.ibzq_qc[iq]
        if self.qsym:
            for i, q in enumerate(self.bzq_qc):
                if abs(q - self.ibzq_qc[iq]).max() < 1e-9:
                    bzq_index = i
                    break    
        else:
            bzq_index = iq
        
        N_c = self.gd.N_c
        eikr_R = np.exp(-2j * pi * np.dot(np.indices(N_c).T, q / N_c).T)

        Gamma = abs(q).max() < 1e-9
        if Gamma and self.skip_gamma:
            return

        Gpk2_G = self.G2_qG[iq]
        if Gamma:
            Gpk2_G = Gpk2_G.copy()
            Gpk2_G[0] = 1.0 / self.gamma

        N = N_c.prod()
        vol = self.gd.dv * N
        nspins = self.nspins

        fcut = 1e-10
        for n1, psit1_R in enumerate(kpt1.psit_nG):
            f1 = kpt1.f_n[n1]
            for n2, psit2_R in enumerate(kpt2.psit_nG):
                if self.acdf:
                    if self.gygi and Gamma:
                        #print n2, kpt2.f_n[n2]/kpt2.weight
                        f2 = (self.q_weights[iq] * kpt2.weight)
                    else:
                        f2 = (self.q_weights[iq] * kpt2.weight
                              * (1 - np.sign(kpt2.eps_n[n2] - kpt1.eps_n[n1])))

                else:
                    f2 = kpt2.f_n[n2] * self.q_weights[iq]
                if abs(f1) < fcut or abs(f2) < fcut:
                    continue
                nt_R = self.calculate_pair_density(n1, n2, kpt1, kpt2,
                                                   ik1, ik2, bzq_index)
                nt_G = self.pwd.fft(nt_R * eikr_R) / N
                vt_G = nt_G.copy()
                vt_G *= -pi * vol / Gpk2_G
                e = np.vdot(nt_G, vt_G).real * nspins * self.hybrid
                self.exx += f1 * f2 * e
                self.exx_kq[ik1,iq] += f1*f2*e        
    
    def calculate_pair_density(self, n1, n2, kpt1, kpt2, ik1, ik2, bzq_index):
        psit1_G = self.kd.transform_wave_function(kpt1.psit_nG[n1], ik1)
        psit2_G = self.kd.transform_wave_function(kpt2.psit_nG[n2], ik2)
        nt_G = psit1_G.conj() * psit2_G
        
        s1 = self.kd.sym_k[ik1]
        s2 = self.kd.sym_k[ik2]
        t1 = self.kd.time_reversal_k[ik1]
        t2 = self.kd.time_reversal_k[ik2]
        k1_c = self.kd.ibzk_kc[kpt1.k]
        k2_c = self.kd.ibzk_kc[kpt2.k]

        Q_aL = {}
        for a in kpt1.P_ani.keys():
            b1 = self.kd.symmetry.a_sa[s1, a]
            b2 = self.kd.symmetry.a_sa[s2, a]
            S1_c = (np.dot(self.spos_ac[a], self.kd.symmetry.op_scc[s1]) -
                   self.spos_ac[b1])
            S2_c = (np.dot(self.spos_ac[a], self.kd.symmetry.op_scc[s2]) -
                   self.spos_ac[b2])
            assert abs(S1_c.round() - S1_c).max() < 1e-13
            assert abs(S2_c.round() - S2_c).max() < 1e-13
            x1 = np.exp(2j * pi * np.dot(k1_c, S1_c))
            x2 = np.exp(2j * pi * np.dot(k2_c, S2_c))
            P1_i = np.dot(self.setups[a].R_sii[s1], kpt1.P_ani[b1][n1]) * x1
            P2_i = np.dot(self.setups[a].R_sii[s2], kpt2.P_ani[b2][n2]) * x2
            if t1:
                P1_i = P1_i.conj()
            if t2:
                P2_i = P2_i.conj()

            D_ii = np.outer(P1_i.conj(), P2_i)
            D_p = pack(D_ii)
            Q_aL[a] = np.dot(D_p, self.setups[a].Delta_pL)

        self.ghat.add(nt_G, Q_aL, bzq_index)
        return nt_G

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
                        p12 = packed_index(i1, i2, ni)
                        exx -= self.hybrid / deg * D_ii[i1, i2] * A

                if setup.X_p is not None:
                    exx -= self.hybrid * np.dot(D_p, setup.X_p)
            exx += self.hybrid * setup.ExxC
        return exx

    def print_initialization(self, xc):
        print('------------------------------------------------------', file=self.txt)
        print('Non-self-consistent HF correlation energy', file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print('Started at:  ', ctime(), file=self.txt)
        print(file=self.txt)
        print('Ground state XC functional     :   %s' % xc, file=self.txt)
        print('Valence electrons              :   %s' % self.setups.nvalence, file=self.txt)
        print('Number of Spins                :   %s' % self.nspins, file=self.txt)
        print('Plane wave cutoff energy       :   %4.1f eV' % (self.ecut*Ha), file=self.txt)
        print('Gamma q-point excluded         :   %s' % self.skip_gamma, file=self.txt)
        if not self.skip_gamma:
            print('Alpha parameter                :   %s' % self.alpha, file=self.txt)
            print('Gamma parameter                :   %3.3f' % self.gamma, file=self.txt)
        print('ACDF method                    :   %s' % self.acdf, file=self.txt)
        print('Number of k-points             :   %s' % len(self.kd.bzk_kc), file=self.txt)
        print('Number of Irreducible k-points :   %s' % len(self.kd.ibzk_kc), file=self.txt)
        print('Number of q-points             :   %s' % len(self.bzq_qc), file=self.txt)
        if not self.qsym:
            print('q-point symmetry               :   %s' % self.qsym, file=self.txt)
        else:
            print('Number of Irreducible q-points :   %s' % len(self.ibzq_qc), file=self.txt)

        print(file=self.txt)
        for q, weight in zip(self.ibzq_qc, self.q_weights):
            print('q: [%1.3f %1.3f %1.3f] - weight: %1.3f' % \
                  (q[0],q[1],q[2], weight/len(self.bzq_qc)), file=self.txt)
        print(file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print('------------------------------------------------------', file=self.txt)
        print(file=self.txt)
        print('Looping over k-points in the full Brillouin zone', file=self.txt)
        print(file=self.txt)
