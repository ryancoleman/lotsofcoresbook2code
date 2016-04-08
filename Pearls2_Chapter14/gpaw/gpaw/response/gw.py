import os
import sys
import pickle
import numpy as np
from math import pi, sqrt
from time import time, ctime
from datetime import timedelta

from ase.parallel import paropen
from ase.units import Hartree, Bohr
from ase.utils import devnull

from gpaw import GPAW
from gpaw.version import version
from gpaw.mpi import world, rank, size, serial_comm
from gpaw.utilities.blas import gemmdot
from gpaw.xc.tools import vxc
from gpaw.wavefunctions.pw import PWWaveFunctions
from gpaw.response.parallel import set_communicator, parallel_partition, SliceAlongFrequency, GatherOrbitals
from gpaw.response.base import BASECHI
from gpaw.response.df0 import DF
from gpaw.response.kernel import calculate_Kxc, calculate_Kc, calculate_Kc_q

class GW(BASECHI):

    def __init__(
                 self,
                 file=None,
                 nbands=None,
                 bands=None,
                 kpoints=None,
                 e_skn=None,
                 eshift=None,
                 w=None,
                 ecut=150.,
                 eta=0.1,
                 ppa=False,
                 E0=None,
                 hilbert_trans=False,
                 wpar=1,
                 vcut=None,
                 numint=False,
                 txt=None
                ):

        BASECHI.__init__(self, calc=file, nbands=nbands, w=w, eshift=eshift, ecut=ecut, eta=eta, txt=devnull)

        self.file = file
        self.gwnbands = nbands
        self.bands = bands
        self.kpoints = kpoints
        self.user_skn = e_skn
        self.ppa = ppa
        self.E0 = E0
        self.hilbert_trans = hilbert_trans
        self.wpar = wpar
        self.vcut = vcut
        self.numint = numint
        self.gwtxtname = txt


    def initialize(self):

        self.ini = True

        BASECHI.initialize(self)

        self.txtname = self.gwtxtname
        self.output_init()

        self.printtxt('GPAW version %s' %(version))
        self.printtxt('-----------------------------------------------')
        self.printtxt('GW calculation started at:')
        self.printtxt(ctime())
        self.printtxt('-----------------------------------------------')
        self.starttime = time()

        calc = self.calc
        kd = self.kd

        # band init
        if self.gwnbands is None:
            if self.npw > calc.wfs.bd.nbands:
                self.nbands = calc.wfs.bd.nbands
            else:
                self.nbands = self.npw

        # eigenvalue init
        if self.user_skn is not None:
            self.printtxt('Use eigenvalues from user.')
            assert np.shape(self.user_skn)[0] == self.nspins, 'Eigenvalues not compatible with .gpw file!'
            assert np.shape(self.user_skn)[1] == self.kd.nibzkpts, 'Eigenvalues not compatible with .gpw file!'
            assert np.shape(self.user_skn)[2] >= self.nbands, 'Too few eigenvalues!'
            self.e_skn = self.user_skn
        else:
            self.printtxt('Use eigenvalues from the calculator.')

        # q point init
        self.bzq_kc = kd.get_bz_q_points()
        self.ibzq_qc = self.bzq_kc # q point symmetry is not used at the moment.
        self.nqpt = np.shape(self.bzq_kc)[0]
        
        # frequency points init
        self.static=False

        if self.ppa: # Plasmon Pole Approximation
            if self.E0 is None:
                self.E0 = Hartree
            self.E0 /= Hartree
            self.w_w = np.array([0., 1j*self.E0])
            self.hilbert_trans=False
            self.wpar=1
        elif self.w_w is None: # static COHSEX
            self.w_w = np.array([0.])
            self.static=True
            self.hilbert_trans=False
            self.wpar=1
            self.eta = 0.0001 / Hartree
        else:
            # create nonlinear frequency grid
            # grid is linear from 0 to wcut with spacing dw
            # spacing is linearily increasing between wcut and wmax
            # Hilbert transforms are still carried out on linear grid
            wcut = self.w_w[0]
            wmax = self.w_w[1]
            dw = self.w_w[2]
            w_w = np.linspace(0., wcut, wcut/dw+1)
            i=1
            wi=wcut
            while wi < wmax:
                wi += i*dw
                w_w = np.append(w_w, wi)
                i+=1
            while len(w_w) % self.wpar != 0:
                wi += i*dw
                w_w = np.append(w_w, wi)
                i+=1

            dw_w = np.zeros(len(w_w))
            dw_w[0] = dw
            dw_w[1:] = w_w[1:] - w_w[:-1]

            self.w_w = w_w / Hartree
            self.dw_w = dw_w / Hartree
            self.eta_w = dw_w * 4 / Hartree
            self.wcut = wcut

            self.wmax = self.w_w[-1]
            self.wmin = self.w_w[0]
            self.dw = self.w_w[1] - self.w_w[0]
            self.Nw = len(self.w_w)
#            self.wpar = int(self.Nw * self.npw**2 * 16. / 1024**2) // 1500 + 1 # estimate memory and parallelize over frequencies

            for s in range(self.nspins):
                emaxdiff = self.e_skn[s][:,self.nbands-1].max() - self.e_skn[s][:,0].min()
                assert (self.wmax > emaxdiff), 'Maximum frequency must be larger than %f' %(emaxdiff*Hartree)

        # GW kpoints init
        if self.kpoints is None:
            self.gwnkpt = self.kd.nibzkpts
            self.gwkpt_k = kd.ibz2bz_k
        else:
            self.gwnkpt = np.shape(self.kpoints)[0]
            self.gwkpt_k = self.kpoints

        # GW bands init
        if self.bands is None:
            self.gwnband = self.nbands
            self.bands = self.gwbands_n = np.arange(self.nbands)
        else:
            self.gwnband = np.shape(self.bands)[0]
            self.gwbands_n = self.bands

        self.alpha = 1j/(2*pi * self.vol * self.kd.nbzkpts)
        
        # parallel init
        assert len(self.w_w) % self.wpar == 0
        self.wcommsize = self.wpar
        self.qcommsize = size // self.wpar
        assert self.qcommsize * self.wcommsize == size, 'wpar must be integer divisor of number of requested cores'
        if self.nqpt != 1: # parallelize over q-points
            self.wcomm, self.qcomm, self.worldcomm = set_communicator(world, rank, size, self.wpar)
            self.ncomm = serial_comm
            self.dfcomm = self.wcomm
            self.kcommsize = 1
        else: # parallelize over bands
            self.wcomm, self.ncomm, self.worldcomm = set_communicator(world, rank, size, self.wpar)
            self.qcomm = serial_comm
            if self.wpar > 1:
                self.dfcomm = self.wcomm
                self.kcommsize = 1
            else:
                self.dfcomm = self.ncomm
                self.kcommsize = self.ncomm.size
        nq, self.nq_local, self.q_start, self.q_end = parallel_partition(
                                  self.nqpt, self.qcomm.rank, self.qcomm.size, reshape=False)
        nb, self.nbands_local, self.m_start, self.m_end = parallel_partition(
                                  self.nbands, self.ncomm.rank, self.ncomm.size, reshape=False)


    def get_QP_spectrum(self, exxfile='EXX.pckl', file='GW.pckl'):

        try:
            self.ini
        except:
            self.initialize()
        self.print_gw_init()
        self.printtxt("calculating Self energy")

        Sigma_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)
        dSigma_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)
        Z_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)

        t0 = time()
        t_w = 0
        t_selfenergy = 0
        for iq in range(self.q_start, self.q_end):
            if iq >= self.nqpt:
                continue
            t1 = time()

            # get screened interaction. 
            df, W_wGG = self.screened_interaction_kernel(iq)
            t2 = time()
            t_w += t2 - t1

            # get self energy
            S, dS = self.get_self_energy(df, W_wGG)
            t3 = time() - t2
            t_selfenergy += t3

            Sigma_skn += S
            dSigma_skn += dS

            del df, W_wGG
            self.timing(iq, t0, self.nq_local, 'iq')

        self.qcomm.barrier()
        self.qcomm.sum(Sigma_skn)
        self.qcomm.sum(dSigma_skn)

        self.printtxt('W_wGG takes %s ' %(timedelta(seconds=round(t_w))))
        self.printtxt('Self energy takes %s ' %(timedelta(seconds=round(t_selfenergy))))

        Z_skn = 1. / (1. - dSigma_skn)

        # exact exchange
        exx = os.path.isfile(exxfile)
        world.barrier()
        if exx:
            open(exxfile)
            self.printtxt("reading Exact exchange and E_XC from file")
        else:
            t0 = time()
            self.get_exact_exchange()
            world.barrier()
            exxfile='EXX.pckl'
            self.printtxt('EXX takes %s ' %(timedelta(seconds=round(time()-t0))))
        data = pickle.load(open(exxfile))
        vxc_skn = data['vxc_skn'] # in Hartree
        exx_skn = data['exx_skn'] # in Hartree
        f_skn = data['f_skn']
        gwkpt_k = data['gwkpt_k']
        gwbands_n = data['gwbands_n']
        assert (gwkpt_k == self.gwkpt_k).all(), 'exxfile inconsistent with input parameters'
        assert (gwbands_n == self.gwbands_n).all(), 'exxfile inconsistent with input parameters'
        if self.user_skn is None:
            e_skn = data['e_skn'] # in Hartree
        else:
            e_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)
        for s in range(self.nspins):
            for i, k in enumerate(self.gwkpt_k):
                ik = self.kd.bz2ibz_k[k]
                for j, n in enumerate(self.gwbands_n):
                    if self.user_skn is None:
                        assert e_skn[s][i][j] == self.e_skn[s][ik][n], 'exxfile inconsistent with eigenvalues'
                    else:
                        e_skn[s][i][j] = self.e_skn[s][ik][n]

        if not self.static:
            QP_skn = e_skn + Z_skn * (Sigma_skn + exx_skn - vxc_skn)
        else:
            QP_skn = e_skn + Sigma_skn - vxc_skn
        self.QP_skn = QP_skn

        # finish
        self.print_gw_finish(e_skn, f_skn, vxc_skn, exx_skn, Sigma_skn, Z_skn, QP_skn)
        data = {
                'gwkpt_k': self.gwkpt_k,
                'gwbands_n': self.gwbands_n,
                'f_skn': f_skn,
                'e_skn': e_skn,         # in Hartree
                'vxc_skn': vxc_skn,     # in Hartree
                'exx_skn': exx_skn,     # in Hartree
                'Sigma_skn': Sigma_skn, # in Hartree
                'Z_skn': Z_skn,         # dimensionless
                'QP_skn': QP_skn        # in Hartree
               }
        if rank == 0:
            pickle.dump(data, open(file, 'w'), -1)


    def screened_interaction_kernel(self, iq):
        """Calcuate W_GG(w) for a given q.
        if static: return W_GG(w=0)
        is not static: return W_GG(q,w) - Vc_GG
        """

        q = self.ibzq_qc[iq]
        w = self.w_w.copy()*Hartree

        optical_limit = False
        if np.abs(q).sum() < 1e-8:
            q = np.array([1e-12, 0, 0]) # arbitrary q, not really need to be calculated
            optical_limit = True

        hilbert_trans = True
        if self.ppa or self.static:
            hilbert_trans = False

        df = DF(calc=self.calc, q=q.copy(), w=w, nbands=self.nbands, eshift=None,
                optical_limit=optical_limit, hilbert_trans=hilbert_trans, xc='RPA', time_ordered=True,
                rpad=self.rpad, vcut=self.vcut, G_plus_q=True,
                eta=self.eta*Hartree, ecut=self.ecut.copy()*Hartree,
                txt='df.out', comm=self.dfcomm, kcommsize=self.kcommsize)

        df.initialize()
        df.e_skn = self.e_skn.copy()
        df.calculate()

        dfinv_wGG = df.get_inverse_dielectric_matrix(xc='RPA')
        assert df.ecut[0] == self.ecut[0]
        if not self.static and not self.ppa:
            assert df.eta == self.eta
            assert df.Nw == self.Nw
            assert df.dw == self.dw

        # calculate Coulomb kernel and use truncation in 2D
        delta_GG = np.eye(df.npw)
        Vq_G = np.diag(calculate_Kc(q,
                                    df.Gvec_Gc,
                                    self.acell_cv,
                                    self.bcell_cv,
                                    self.pbc,
                                    integrate_gamma=True,
                                    N_k=self.kd.N_c,
                                    vcut=self.vcut))**0.5
        if (self.vcut == '2D' and df.optical_limit) or self.numint:
            for iG in range(len(df.Gvec_Gc)):
                if df.Gvec_Gc[iG, 0] == 0 and df.Gvec_Gc[iG, 1] == 0:
                    v_q, v0_q = calculate_Kc_q(self.acell_cv,
                                               self.bcell_cv,
                                               self.pbc,
                                               self.kd.N_c,
                                               vcut=self.vcut,
                                               q_qc=np.array([q]),
                                               Gvec_c=df.Gvec_Gc[iG])
                    Vq_G[iG] = v_q[0]**0.5
        self.Kc_GG = np.outer(Vq_G, Vq_G)

        if self.ppa:
            dfinv1_GG = dfinv_wGG[0] - delta_GG
            dfinv2_GG = dfinv_wGG[1] - delta_GG
            self.wt_GG = self.E0 * np.sqrt(dfinv2_GG / (dfinv1_GG - dfinv2_GG))
            self.R_GG = - self.wt_GG / 2 * dfinv1_GG
            del dfinv_wGG
            dfinv_wGG = np.array([1j*pi*self.R_GG + delta_GG])

        if self.static:
            assert len(dfinv_wGG) == 1
            W_GG = dfinv_wGG[0] * self.Kc_GG

            return df, W_GG
        else:
            Nw = np.shape(dfinv_wGG)[0]
            W_wGG = np.zeros_like(dfinv_wGG)
            for iw in range(Nw):
                dfinv_wGG[iw] -= delta_GG
                W_wGG[iw] = dfinv_wGG[iw] * self.Kc_GG

            return df, W_wGG


    def get_self_energy(self, df, W_wGG):

        Sigma_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)
        dSigma_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)

        wcomm = df.wcomm

        if self.static:
            W_wGG = np.array([W_wGG])

        if not self.hilbert_trans: #method 1
            Wbackup_wG0 = W_wGG[:,:,0].copy()
            Wbackup_w0G = W_wGG[:,0,:].copy()

        else: #method 2, perform Hilbert transform
            nG = np.shape(W_wGG)[1]
            coords = np.zeros(wcomm.size, dtype=int)
            nG_local = nG**2 // wcomm.size
            if wcomm.rank == wcomm.size - 1:
                  nG_local = nG**2 - (wcomm.size - 1) * nG_local
            wcomm.all_gather(np.array([nG_local]), coords)
            W_Wg = SliceAlongFrequency(W_wGG, coords, wcomm)

            ng = np.shape(W_Wg)[1]
            Nw = int(self.w_w[-1] / self.dw)

            w1_ww = np.zeros((Nw, df.Nw), dtype=complex)
            for iw in range(Nw):
                w1 = iw * self.dw
                w1_ww[iw] = 1./(w1 + self.w_w + 1j*self.eta_w) + 1./(w1 - self.w_w + 1j*self.eta_w)
                w1_ww[iw,0] -= 1./(w1 + 1j*self.eta_w[0]) # correct w'=0
                w1_ww[iw] *= self.dw_w

            Cplus_Wg = np.zeros((Nw, ng), dtype=complex)
            Cminus_Wg = np.zeros((Nw, ng), dtype=complex)
            Cplus_Wg = gemmdot(w1_ww, W_Wg, beta=0.0)
            Cminus_Wg = gemmdot(w1_ww.conj(), W_Wg, beta=0.0)

        for s in range(self.nspins):
            for i, k in enumerate(self.gwkpt_k): # k is bzk index
                if df.optical_limit:
                    kq_c = df.kd.bzk_kc[k]
                else:
                    kq_c = df.kd.bzk_kc[k] - df.q_c  # k - q

                kq = df.kd.where_is_q(kq_c, df.kd.bzk_kc)            
                assert df.kq_k[kq] == k
                ibzkpt1 = df.kd.bz2ibz_k[k]
                ibzkpt2 = df.kd.bz2ibz_k[kq]

                for j, n in enumerate(self.bands):
                    for m in range(self.m_start, self.m_end):

                        if self.e_skn[s][ibzkpt2, m] > self.eFermi:
                            sign = 1.
                        else:
                            sign = -1.

                        rho_G = df.density_matrix(m, n, kq, spin1=s, spin2=s)

                        if not self.hilbert_trans: #method 1
                            W_wGG[:,:,0] = Wbackup_wG0
                            W_wGG[:,0,:] = Wbackup_w0G

                            # w1 = w - epsilon_m,k-q
                            w1 = self.e_skn[s][ibzkpt1,n] - self.e_skn[s][ibzkpt2,m]

                            if self.ppa:
                                # analytical expression for Plasmon Pole Approximation
                                W_GG = sign * W_wGG[0] * (1./(w1 + self.wt_GG - 1j*self.eta) -
                                                          1./(w1 - self.wt_GG + 1j*self.eta))
                                W_GG -= W_wGG[0] * (1./(w1 + self.wt_GG + 1j*self.eta*sign) +
                                                    1./(w1 - self.wt_GG + 1j*self.eta*sign))
                                W_G = gemmdot(W_GG, rho_G, beta=0.0)
                                Sigma_skn[s,i,j] += np.real(gemmdot(W_G, rho_G, alpha=self.alpha, beta=0.0,trans='c'))

                                W_GG = sign * W_wGG[0] * (1./(w1 - self.wt_GG + 1j*self.eta)**2 -
                                                          1./(w1 + self.wt_GG - 1j*self.eta)**2)
                                W_GG += W_wGG[0] * (1./(w1 - self.wt_GG + 1j*self.eta*sign)**2 +
                                                    1./(w1 + self.wt_GG + 1j*self.eta*sign)**2)
                                W_G = gemmdot(W_GG, rho_G, beta=0.0)
                                dSigma_skn[s,i,j] += np.real(gemmdot(W_G, rho_G, alpha=self.alpha, beta=0.0,trans='c'))

                            elif self.static:
                                W1_GG = W_wGG[0] - np.eye(df.npw)*self.Kc_GG
                                W2_GG = W_wGG[0]

                                # perform W_GG * np.outer(rho_G.conj(), rho_G).sum(GG)
                                W_G = gemmdot(W1_GG, rho_G, beta=0.0) # Coulomb Hole
                                Sigma_skn[s,i,j] += np.real(gemmdot(W_G, rho_G, alpha=self.alpha*pi/1j, beta=0.0,trans='c'))
                                if sign == -1:
                                    W_G = gemmdot(W2_GG, rho_G, beta=0.0) # Screened Exchange
                                    Sigma_skn[s,i,j] -= np.real(gemmdot(W_G, rho_G, alpha=2*self.alpha*pi/1j, beta=0.0,trans='c'))
                                del W1_GG, W2_GG, W_G, rho_G

                            else:
                                # perform W_wGG * np.outer(rho_G.conj(), rho_G).sum(GG)
                                W_wG = gemmdot(W_wGG, rho_G, beta=0.0)
                                C_wlocal = gemmdot(W_wG, rho_G, alpha=self.alpha, beta=0.0,trans='c')
                                del W_wG, rho_G

                                C_w = np.zeros(df.Nw, dtype=complex)
                                wcomm.all_gather(C_wlocal, C_w)
                                del C_wlocal

                                # calculate self energy
                                w1_w = 1./(w1 - self.w_w + 1j*self.eta_w*sign) + 1./(w1 + self.w_w + 1j*self.eta_w*sign)
                                w1_w[0] -= 1./(w1 + 1j*self.eta_w[0]*sign) # correct w'=0
                                w1_w *= self.dw_w
                                Sigma_skn[s,i,j] += np.real(gemmdot(C_w, w1_w, beta=0.0))

                                # calculate derivate of self energy with respect to w
                                w1_w = 1./(w1 - self.w_w + 1j*self.eta_w*sign)**2 + 1./(w1 + self.w_w + 1j*self.eta_w*sign)**2
                                w1_w[0] -= 1./(w1 + 1j*self.eta_w[0]*sign)**2 # correct w'=0
                                w1_w *= self.dw_w
                                dSigma_skn[s,i,j] -= np.real(gemmdot(C_w, w1_w, beta=0.0))

                        else: #method 2
                            if not np.abs(self.e_skn[s][ibzkpt2,m] - self.e_skn[s][ibzkpt1,n]) < 1e-10:
                                sign *= np.sign(self.e_skn[s][ibzkpt1,n] - self.e_skn[s][ibzkpt2,m])

                            # find points on frequency grid
                            w0 = self.e_skn[s][ibzkpt1,n] - self.e_skn[s][ibzkpt2,m]
                            w0_id = np.abs(int(w0 / self.dw))
                            w1 = w0_id * self.dw
                            w2 = (w0_id + 1) * self.dw

                            # choose plus or minus, treat optical limit:
                            if sign == 1:
                                C_Wg = Cplus_Wg[w0_id:w0_id+2] # only two grid points needed for each w0
                            if sign == -1:
                                C_Wg = Cminus_Wg[w0_id:w0_id+2] # only two grid points needed for each w0

                            C_wGG = GatherOrbitals(C_Wg, coords, wcomm).copy()
                            del C_Wg

                            # special treat of w0 = 0 (degenerate states):
                            if w0_id == 0:
                                Cplustmp_GG = GatherOrbitals(Cplus_Wg[1], coords, wcomm).copy()
                                Cminustmp_GG = GatherOrbitals(Cminus_Wg[1], coords, wcomm).copy()

                            # perform C_wGG * np.outer(rho_G.conj(), rho_G).sum(GG)

                            if w0_id == 0:
                                Sw0_G = gemmdot(C_wGG[0], rho_G, beta=0.0)
                                Sw0 = np.real(gemmdot(Sw0_G, rho_G, alpha=self.alpha, beta=0.0, trans='c'))
                                Sw1_G = gemmdot(Cplustmp_GG, rho_G, beta=0.0)
                                Sw1 = np.real(gemmdot(Sw1_G, rho_G, alpha=self.alpha, beta=0.0, trans='c'))
                                Sw2_G = gemmdot(Cminustmp_GG, rho_G, beta=0.0)
                                Sw2 = np.real(gemmdot(Sw2_G, rho_G, alpha=self.alpha, beta=0.0, trans='c'))

                                Sigma_skn[s,i,j] += Sw0
                                dSigma_skn[s,i,j] += (Sw1 + Sw2)/(2*self.dw)

                            else:                        
                                Sw1_G = gemmdot(C_wGG[0], rho_G, beta=0.0)
                                Sw1 = np.real(gemmdot(Sw1_G, rho_G, alpha=self.alpha, beta=0.0, trans='c'))
                                Sw2_G = gemmdot(C_wGG[1], rho_G, beta=0.0)
                                Sw2 = np.real(gemmdot(Sw2_G, rho_G, alpha=self.alpha, beta=0.0, trans='c'))

                                Sw0 = (w2-np.abs(w0))/self.dw * Sw1 + (np.abs(w0)-w1)/self.dw * Sw2
                                Sigma_skn[s,i,j] += np.sign(self.e_skn[s][ibzkpt1,n] - self.e_skn[s][ibzkpt2,m]) * Sw0
                                dSigma_skn[s,i,j] += (Sw2 - Sw1)/self.dw

        self.ncomm.barrier()
        self.ncomm.sum(Sigma_skn)
        self.ncomm.sum(dSigma_skn)

        return Sigma_skn, dSigma_skn 


    def get_exact_exchange(self, ecut=None, communicator=world, file='EXX.pckl'):

        try:
            self.ini
        except:
            self.initialize()

        self.printtxt('------------------------------------------------')
        self.printtxt("calculating Exact exchange and E_XC")

        calc = GPAW(self.file, communicator=communicator, parallel={'domain':1, 'band':1}, txt=None)
        v_xc = vxc(calc)

        if ecut is None:
            ecut = self.ecut.max()
        else:
            ecut /= Hartree

        if self.pwmode: # planewave mode
            from gpaw.xc.hybridg import HybridXC
            self.printtxt('Use planewave ecut from groundstate calculator: %4.1f eV' % (calc.wfs.pd.ecut*Hartree) )
            exx = HybridXC('EXX', alpha=5.0, bandstructure=True, bands=self.bands)
        else:                                     # grid mode
            from gpaw.xc.hybridk import HybridXC
            self.printtxt('Planewave ecut (eV): %4.1f' % (ecut*Hartree) )
            exx = HybridXC('EXX', alpha=5.0, ecut=ecut, bands=self.bands)
        calc.get_xc_difference(exx)

        e_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)
        f_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)
        vxc_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)
        exx_skn = np.zeros((self.nspins, self.gwnkpt, self.gwnband), dtype=float)

        for s in range(self.nspins):
            for i, k in enumerate(self.gwkpt_k):
                ik = self.kd.bz2ibz_k[k]
                for j, n in enumerate(self.gwbands_n):
                    e_skn[s][i][j] = self.e_skn[s][ik][n]
                    f_skn[s][i][j] = self.f_skn[s][ik][n]
                    vxc_skn[s][i][j] = v_xc[s][ik][n] / Hartree
                    exx_skn[s][i][j] = exx.exx_skn[s][ik][n]
                    if self.eshift is not None:
                        if e_skn[s][i][j] > self.eFermi:
                            vxc_skn[s][i][j] += self.eshift / Hartree

        data = {
                'e_skn':     e_skn,      # in Hartree
                'vxc_skn':   vxc_skn,    # in Hartree
                'exx_skn':   exx_skn,    # in Hartree
                'f_skn':     f_skn,
                'gwkpt_k':   self.gwkpt_k,
                'gwbands_n': self.gwbands_n
               }
        if rank == 0:
            pickle.dump(data, open(file, 'w'), -1)
            self.printtxt("------------------------------------------------")
            self.printtxt("non-selfconsistent HF eigenvalues are (eV):")
            self.printtxt((e_skn - vxc_skn + exx_skn)*Hartree)


    def print_gw_init(self):

        if self.eshift is not None:
            self.printtxt('Shift unoccupied bands by %f eV' % (self.eshift))
        for s in range(self.nspins):
            self.printtxt('Lowest eigenvalue (spin=%s) : %f eV' %(s, self.e_skn[s][:, 0].min()*Hartree))
            self.printtxt('Highest eigenvalue (spin=%s): %f eV' %(s, self.e_skn[s][:, self.nbands-1].max()*Hartree))
        self.printtxt('')
        if self.ecut[0] == self.ecut[1] and self.ecut[0] == self.ecut[2]:
            self.printtxt('Plane wave ecut (eV)         : %4.1f' % (self.ecut[0]*Hartree))
        else:
            self.printtxt('Plane wave ecut (eV)         : (%f, %f, %f)' % (self.ecut[0]*Hartree,self.ecut[1]*Hartree,self.ecut[2]*Hartree) )
        self.printtxt('Number of plane waves used   : %d' %(self.npw) )
        self.printtxt('Number of bands              : %d' %(self.nbands) )
        self.printtxt('Number of k points           : %d' %(self.kd.nbzkpts) )
        self.printtxt("Number of IBZ k points       : %d" %(self.kd.nibzkpts))
        self.printtxt("Number of spins              : %d" %(self.nspins))
        self.printtxt('')
        if self.ppa:
            self.printtxt("Use Plasmon Pole Approximation")
            self.printtxt("imaginary frequency (eV)     : %.2f" %(self.E0*Hartree))
            self.printtxt("broadening (eV)              : %.2f" %(self.eta*Hartree))
        elif self.static:
            self.printtxt("Use static COHSEX")
        else:
            self.printtxt("Linear frequency grid (eV)   : %.2f - %.2f in %.2f" %(self.wmin*Hartree, self.wcut, self.dw*Hartree))
            self.printtxt("Maximum frequency (eV)       : %.2f" %(self.wmax*Hartree))
            self.printtxt("Number of frequency points   : %d" %(self.Nw))
            self.printtxt("Use Hilbert transform        : %s" %(self.hilbert_trans))
        self.printtxt('')
        self.printtxt('Coulomb interaction cutoff   : %s' % self.vcut)
        self.printtxt('')
        self.printtxt('Calculate matrix elements for k = :')
        for k in self.gwkpt_k:
            self.printtxt(self.kd.bzk_kc[k])
        self.printtxt('')
        self.printtxt('Calculate matrix elements for n = :')
        self.printtxt(self.gwbands_n)
        self.printtxt('')


    def print_gw_finish(self, e_skn, f_skn, vxc_skn, exx_skn, Sigma_skn, Z_skn, QP_skn):

        self.printtxt("------------------------------------------------")
        self.printtxt("Kohn-Sham eigenvalues are (eV): ")
        self.printtxt("%s \n" %(e_skn*Hartree))
        self.printtxt("Occupation numbers are: ")
        self.printtxt("%s \n" %(f_skn*self.kd.nbzkpts))
        self.printtxt("Kohn-Sham exchange-correlation contributions are (eV): ")
        self.printtxt("%s \n" %(vxc_skn*Hartree))
        self.printtxt("Exact exchange contributions are (eV): ")
        self.printtxt("%s \n" %(exx_skn*Hartree))
        self.printtxt("Self energy contributions are (eV):")
        self.printtxt("%s \n" %(Sigma_skn*Hartree))
        if not self.static:
            self.printtxt("Renormalization factors are:")
            self.printtxt("%s \n" %(Z_skn))

        totaltime = round(time() - self.starttime)
        self.printtxt("GW calculation finished in %s " %(timedelta(seconds=totaltime)))
        self.printtxt("------------------------------------------------")
        self.printtxt("Quasi-particle energies are (eV): ")
        self.printtxt(QP_skn*Hartree)
