from __future__ import print_function
import sys
from time import time, ctime
import numpy as np
from math import sqrt, pi
from datetime import timedelta

from ase.units import Hartree, Bohr
from ase.utils import devnull

from gpaw import GPAW, extra_parameters
from gpaw.utilities import unpack
from gpaw.utilities.blas import gemmdot, gemv
from gpaw.mpi import world, rank, size, serial_comm
from gpaw.lfc import LocalizedFunctionsCollection as LFC
from gpaw.grid_descriptor import GridDescriptor
from gpaw.utilities.memory import maxrss
from gpaw.fd_operators import Gradient
from gpaw.response.cell import get_primitive_cell, set_Gvectors
from gpaw.response.math_func import delta_function,  \
     two_phi_planewave_integrals
from gpaw.response.parallel import set_communicator, \
     parallel_partition, SliceAlongFrequency, SliceAlongOrbitals
from gpaw.response.kernel import calculate_Kxc, calculate_Kc, calculate_Kc_q
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.wavefunctions.pw import PWLFC
import gpaw.wavefunctions.pw as pw

class BASECHI:
    """This class is to store the basic common stuff for chi and bse."""

    def __init__(self,
                 calc=None,
                 nbands=None,
                 w=None,
                 q=None,
                 eshift=None,
                 ecut=10.,
                 density_cut=None,
                 G_plus_q=False,
                 eta=0.2,
                 rpad=None,
                 ftol=1e-5,
                 txt=None,
                 optical_limit=False):

        if rpad is None:
            rpad = np.ones(3, int)

        self.txtname = txt
        self.output_init()

        if isinstance(calc, str):
            # Always use serial_communicator when a filename is given.
            self.calc = GPAW(calc, communicator=serial_comm, txt=None)
        else:
            # To be optimized so that the communicator is loaded automatically 
            # according to kcommsize.
            # 
            # so temporarily it is used like this :
            # kcommsize = int (should <= world.size)
            # r0 = rank % kcommsize
            # ranks = np.arange(r0, r0+size, kcommsize)
            # calc = GPAW(filename.gpw, communicator=ranks, txt=None)
            self.calc = calc
        
        if self.calc is not None:
            self.pwmode = isinstance(self.calc.wfs, pw.PWWaveFunctions)
        else:
            self.pwmode = False
        if self.pwmode:
            assert self.calc.wfs.world.size == 1

        self.nbands = nbands
        self.q_c = q

        # chi.py modifies the input array w by dividing by Hartree.
        # This will change the user-supplied arrays in-place unless
        # we create a copy.  So now we create a copy.  *Grumble*
        #
        # To make matters worse, w is allowed to be None (why not take
        # care of that *before*??  This should really be cleaned up.
        if isinstance(w, np.ndarray):
            w = w.copy()
        self.w_w = w
        self.eta = eta
        self.ftol = ftol
        if isinstance(ecut, int) or isinstance(ecut, float):
            self.ecut = np.ones(3) * ecut
        else:
            assert len(ecut) == 3
            self.ecut = np.array(ecut, dtype=float)
        self.density_cut = density_cut
        self.G_plus_q = G_plus_q
        self.rpad = rpad
        self.optical_limit = optical_limit
        if self.optical_limit:
            self.qopt = 1e-5
        self.eshift = eshift

    def initialize(self):
                        
        self.eta /= Hartree
        self.ecut /= Hartree

        calc = self.calc
        self.nspins = self.calc.wfs.nspins

        # kpoint init
        self.kd = kd = calc.wfs.kd
        self.nikpt = kd.nibzkpts
        self.ftol /= kd.nbzkpts

        # cell init
        self.acell_cv = calc.wfs.gd.cell_cv
        self.acell_cv, self.bcell_cv, self.vol, self.BZvol = \
                       get_primitive_cell(self.acell_cv,rpad=self.rpad)

        # grid init
        gd = calc.wfs.gd.new_descriptor(comm=serial_comm)
        self.pbc = gd.pbc_c
        self.gd = gd
        self.nG0 = np.prod(gd.N_c)
        # Number of grid points and volume including zero padding
        self.nGrpad = gd.N_c * self.rpad
        self.nG0rpad = np.prod(self.nGrpad)
        self.d_c = [Gradient(gd, i, n=4, dtype=complex).apply for i in range(3)]

        # obtain eigenvalues, occupations
        nibzkpt = kd.nibzkpts
        kweight_k = kd.weight_k

        self.eFermi = self.calc.occupations.get_fermi_level()

        try:
            self.e_skn
            self.printtxt('Use eigenvalues from user.')
        except:
            self.printtxt('Use eigenvalues from the calculator.')
            self.e_skn = {}
            self.f_skn = {}
            for ispin in range(self.nspins):
                self.e_skn[ispin] = np.array([calc.get_eigenvalues(kpt=k, spin=ispin)
                                              for k in range(nibzkpt)]) / Hartree
                self.f_skn[ispin] = np.array([calc.get_occupation_numbers(kpt=k, spin=ispin)
                                              / kweight_k[k]
                                              for k in range(nibzkpt)]) / kd.nbzkpts
            #self.printtxt('Eigenvalues(k=0) are:')
            #print  >> self.txt, self.e_skn[0][0] * Hartree

        self.enoshift_skn = {}
        for ispin in range(self.nspins):
            self.enoshift_skn[ispin] = self.e_skn[ispin].copy()
        if self.eshift is not None:
            self.add_discontinuity(self.eshift)
            self.printtxt('Shift unoccupied bands by %f eV' % (self.eshift))
        # k + q init
        if self.q_c is not None:
            self.qq_v = np.dot(self.q_c, self.bcell_cv) # summation over c
    
            if self.optical_limit:
                kq_k = np.arange(kd.nbzkpts)
                self.expqr_g = 1.
            else:
                r_vg = gd.get_grid_point_coordinates() # (3, nG)
                qr_g = gemmdot(self.qq_v, r_vg, beta=0.0)
                self.expqr_g = np.exp(-1j * qr_g)
                del r_vg, qr_g
                kq_k = kd.find_k_plus_q(self.q_c)
            self.kq_k = kq_k

        # Plane wave init
        if self.G_plus_q:
            self.npw, self.Gvec_Gc, self.Gindex_G = set_Gvectors(self.acell_cv,
                                                                 self.bcell_cv,
                                                                 self.gd.N_c,
                                                                 self.ecut,
                                                                 q=self.q_c)
        else:
            self.npw, self.Gvec_Gc, self.Gindex_G = set_Gvectors(self.acell_cv,
                                                                 self.bcell_cv,
                                                                 self.gd.N_c,
                                                                 self.ecut)

        # band init
        if self.nbands is None:
            self.nbands = calc.wfs.bd.nbands
        self.nvalence = calc.wfs.nvalence

        # Projectors init
        setups = calc.wfs.setups
        self.spos_ac = calc.atoms.get_scaled_positions()

        if self.pwmode:
            self.pt = PWLFC([setup.pt_j for setup in setups], self.calc.wfs.pd)
            self.pt.set_positions(self.spos_ac)
        else:
            self.pt = LFC(gd, [setup.pt_j for setup in setups],
                          KPointDescriptor(self.kd.bzk_kc),
                          dtype=complex, forces=True)

            self.pt.set_positions(self.spos_ac)

        # Printing calculation information
        self.print_stuff()

        return

    def output_init(self):

        if self.txtname is None:
            if rank == 0:
                self.txt = sys.stdout
            else:
                sys.stdout = devnull
                self.txt = devnull
        elif self.txtname == devnull:
            self.txt = devnull
        else:
            assert type(self.txtname) is str
            from ase.parallel import paropen
            self.txt = paropen(self.txtname,'w')


    def printtxt(self, text):
        print(text, file=self.txt)


    def print_stuff(self):

        printtxt = self.printtxt
        printtxt('')
        printtxt('Parameters used:')
        printtxt('')
        printtxt('Unit cell (a.u.):')
        printtxt(self.acell_cv)
        printtxt('Volume of cell (a.u.**3)     : %f' % self.vol)
        printtxt('Reciprocal cell (1/a.u.)')
        printtxt(self.bcell_cv)
        printtxt('BZ volume (1/a.u.**3)        : %f' % self.BZvol)
        printtxt('Number of G-vectors / Grid   : %d %s'
                 % (self.nG0, tuple(self.gd.N_c)))
        printtxt('')                         
        printtxt('Coulomb interaction cutoff   : %s' % self.vcut)                            
        printtxt('')                         
        printtxt('Number of bands              : %d' % self.nbands)
        printtxt('Number of kpoints            : %d' % self.kd.nbzkpts)
        if self.ecut[0] == self.ecut[1] and self.ecut[0] == self.ecut[2]:
            printtxt('Planewave ecut (eV)          : %4.1f' % (self.ecut[0] * Hartree))
        else:
            printtxt('Planewave ecut (eV)          : (%f, %f, %f)' % tuple(self.ecut * Hartree))
        printtxt('Number of planewave used     : %d' % self.npw)
        printtxt('Broadening (eta)             : %f' % (self.eta * Hartree))
        printtxt('')
        if self.q_c is not None:
            if self.optical_limit:
                printtxt('Optical limit calculation ! (q=1e-5)')
            else:
                printtxt('q in reduced coordinate        : (%f %f %f)' % tuple(self.q_c))
                printtxt('q in cartesian coordinate (1/A): (%f %f %f)' % tuple(self.qq_v / Bohr))
                printtxt('|q| (1/A)                      : %f' % np.linalg.norm(self.qq_v / Bohr))


    def timing(self, i, t0, n_local, txt):
        if i == 0:
            dt = time() - t0
            self.totaltime = dt * n_local
            self.printtxt('  Finished %s 0 in %s, estimate %s left.'
                          % (txt, timedelta(seconds=round(dt)),
                             timedelta(seconds=round(self.totaltime))))
        if rank == 0 and n_local // 5 > 0:            
            if i > 0 and i % (n_local // 5) == 0:
                dt =  time() - t0
                self.printtxt('  Finished %s %d in %s, estimate %s left.'
                              % (txt, i, timedelta(seconds=round(dt)),
                                 timedelta(seconds=round(self.totaltime
                                                         - dt))))

    def get_phi_aGp(self, q_c=None, parallel=True, alldir=False):
        if q_c is None:
            q_c = self.q_c
            qq_v = self.qq_v
            optical_limit = self.optical_limit
        else:
            optical_limit = False
            if np.abs(q_c).sum() < 1e-8:
                q_c = np.array([0.0001, 0, 0])
                optical_limit = True
            qq_v = np.dot(q_c, self.bcell_cv)
            
        setups = self.calc.wfs.setups
        spos_ac = self.calc.atoms.get_scaled_positions()
        
        kk_Gv = gemmdot(q_c + self.Gvec_Gc, self.bcell_cv.copy(), beta=0.0)
        phi_aGp = {}
        phiG0_avp = {}

        if parallel:
            from gpaw.response.parallel import parallel_partition

            npw, npw_local, Gstart, Gend = parallel_partition(
                               self.npw, self.comm.rank, self.comm.size, reshape=False)
        else:
            Gstart = 0
            Gend = self.npw
        
        for a, id in enumerate(setups.id_a):
            phi_aGp[a] = two_phi_planewave_integrals(kk_Gv, setups[a], Gstart, Gend)
            for iG in range(Gstart, Gend):
                phi_aGp[a][iG] *= np.exp(-1j * 2. * pi *
                                         np.dot(q_c + self.Gvec_Gc[iG], spos_ac[a]) )
            if parallel:
                self.comm.sum(phi_aGp[a])
        # For optical limit, G == 0 part should change
        if optical_limit:
            for a, id in enumerate(setups.id_a):
                nabla_iiv = setups[a].nabla_iiv
                phi_aGp[a][0] = -1j * (np.dot(nabla_iiv, qq_v)).ravel()

                phiG0_avp[a] = np.zeros((3, len(phi_aGp[a][0])), complex)
                for dir in range(3): # 3 dimension
                    q2_c = np.diag((1,1,1))[dir] * self.qopt
                    qq2_v = np.dot(q2_c, self.bcell_cv) # summation over c
                    phiG0_avp[a][dir] = -1j * (np.dot(nabla_iiv, qq2_v)).ravel()

        if alldir:
            return phi_aGp, phiG0_avp
        else:
            return phi_aGp


    def get_wavefunction(self, ibzk, n, check_focc=True, spin=0):
        if (self.calc.wfs.world.size == 1 or self.calc.wfs.gd.comm.size != 1 
            or self.calc.input_parameters['mode'] == 'lcao'):
            if not check_focc:
                return
            else:
                psit_G = self.calc.wfs.get_wave_function_array(n, ibzk, spin)
            
                if self.calc.wfs.world.size == 1:
                    return np.complex128(psit_G)
                
                if self.calc.wfs.world.rank != 0:
                    psit_G = self.calc.wfs.gd.empty(dtype=self.calc.wfs.dtype,
                                                    global_array=True)
                self.calc.wfs.world.broadcast(psit_G, 0)
        
                return np.complex128(psit_G)
        else:
            # support ground state calculation with kpoint and band parallelization
            # but domain decomposition must = 1
            kpt_rank, u = self.calc.wfs.kd.get_rank_and_index(0, ibzk)
            bzkpt_rank = self.kcomm.rank
            band_rank, myn = self.calc.wfs.bd.who_has(n)
            assert self.calc.wfs.gd.comm.size == 1
            world_rank = (kpt_rank * self.calc.wfs.band_comm.size + band_rank)

            # in the following, kpt_rank is assigned to world_rank
            klist = np.array([world_rank, u, bzkpt_rank, myn])
            klist_kcomm = np.zeros((self.kcomm.size, 4), dtype=int)            
            self.kcomm.all_gather(klist, klist_kcomm)

            check_focc_global = np.zeros(self.kcomm.size, dtype=bool)
            self.kcomm.all_gather(np.array([check_focc]), check_focc_global)

            psit_G = self.calc.wfs.gd.empty(dtype=self.calc.wfs.dtype)

            for i in range(self.kcomm.size):
                if check_focc_global[i]:
                    kpt_rank, u, bzkpt_rank, nlocal = klist_kcomm[i]
                    if kpt_rank == bzkpt_rank:
                        if rank == kpt_rank:
                            psit_G = self.calc.wfs.kpt_u[u].psit_nG[nlocal]
                    else:
                        if rank == kpt_rank:
                            world.send(self.calc.wfs.kpt_u[u].psit_nG[nlocal],
                                       bzkpt_rank, 1300+bzkpt_rank)
                        if rank == bzkpt_rank:
                            psit_G = self.calc.wfs.gd.empty(dtype=self.calc.wfs.dtype)
                            world.receive(psit_G, kpt_rank, 1300+bzkpt_rank)
                    
            self.wScomm.broadcast(psit_G, 0)

            return psit_G


    def add_discontinuity(self, shift):

        for ispin in range(self.nspins):
            for k in range(self.kd.nibzkpts):
                for i in range(self.e_skn[0].shape[1]):
                    if self.e_skn[ispin][k,i] > self.eFermi:
                        self.e_skn[ispin][k,i] += shift / Hartree

    def density_matrix(self, n, m, k, kq=None,
                       spin1=0, spin2=0, phi_aGp=None, Gspace=True):

        gd = self.gd
        kd = self.kd
        optical_limit = False

        if kq is None:
            kq = self.kq_k[k]
            expqr_g = self.expqr_g
            q_v = self.qq_v
            optical_limit = self.optical_limit
            q_c = self.q_c
        else:
            q_c = kd.bzk_kc[kq] - kd.bzk_kc[k]
            q_c[np.where(q_c>0.501)] -= 1
            q_c[np.where(q_c<-0.499)] += 1
            
            if (np.abs(q_c) < self.ftol).all():
                optical_limit = True
                q_c = self.q_c
            q_v = np.dot(q_c, self.bcell_cv)
            r_vg = gd.get_grid_point_coordinates() # (3, nG)
            qr_g = gemmdot(q_v, r_vg, beta=0.0)
            expqr_g = np.exp(-1j * qr_g)
            if optical_limit:
                expqr_g = 1

        ibzkpt1 = kd.bz2ibz_k[k]
        ibzkpt2 = kd.bz2ibz_k[kq]
        
        psitold_g = self.get_wavefunction(ibzkpt1, n, True, spin=spin1)
        psit1_g = kd.transform_wave_function(psitold_g, k)
        
        psitold_g = self.get_wavefunction(ibzkpt2, m, True, spin=spin2)
        psit2_g = kd.transform_wave_function(psitold_g, kq)

        if Gspace is False:
            return psit1_g.conj() * psit2_g * expqr_g
        else:
            tmp_g = psit1_g.conj()* psit2_g * expqr_g
            # zero padding is included through the FFT
            rho_g = np.fft.fftn(tmp_g, s=self.nGrpad) * self.vol / self.nG0rpad

            # Here, planewave cutoff is applied
            rho_G = rho_g.ravel()[self.Gindex_G]

            if optical_limit:
                dpsit_g = gd.empty(dtype=complex)
                tmp = np.zeros((3), dtype=complex)
                phase_cd = np.exp(2j * pi * gd.sdisp_cd * kd.bzk_kc[kq, :, np.newaxis])
                for ix in range(3):
                    self.d_c[ix](psit2_g, dpsit_g, phase_cd)
                    tmp[ix] = gd.integrate(psit1_g.conj() * dpsit_g)
                rho_G[0] = -1j * np.dot(q_v, tmp)
                
            calc = self.calc
            pt = self.pt
            if not self.pwmode:
                if calc.wfs.world.size > 1 or kd.nbzkpts == 1:
                    P1_ai = pt.dict()
                    pt.integrate(psit1_g, P1_ai, k)
                    P2_ai = pt.dict()
                    pt.integrate(psit2_g, P2_ai, kq)
                else:
                    P1_ai = self.get_P_ai(k, n, spin1)
                    P2_ai = self.get_P_ai(kq, m, spin2)
            else:
                # first calculate P_ai at ibzkpt, then rotate to k
                u = self.kd.get_rank_and_index(spin1, ibzkpt1)[1]
                Ptmp_ai = pt.dict()
                kpt = calc.wfs.kpt_u[u]
                pt.integrate(kpt.psit_nG[n], Ptmp_ai, ibzkpt1)
                P1_ai = self.get_P_ai(k, n, spin1, Ptmp_ai)

                u = self.kd.get_rank_and_index(spin2, ibzkpt2)[1]
                Ptmp_ai = pt.dict()
                kpt = calc.wfs.kpt_u[u]
                pt.integrate(kpt.psit_nG[m], Ptmp_ai, ibzkpt2)
                P2_ai = self.get_P_ai(kq, m, spin2, Ptmp_ai)

            if phi_aGp is None:
                try:
                    if not self.mode == 'RPA':
                        if optical_limit:
                            iq = kd.where_is_q(np.zeros(3), self.bzq_qc)
                        else:
                            iq = kd.where_is_q(q_c, self.bzq_qc)
                            assert np.abs(self.bzq_qc[iq] - q_c).sum() < 1e-8

                    phi_aGp = self.load_phi_aGp(self.reader, iq) #phi_qaGp[iq]
                except AttributeError:
                    phi_aGp = self.phi_aGp

            for a, id in enumerate(self.calc.wfs.setups.id_a):
                P_p = np.outer(P1_ai[a].conj(), P2_ai[a]).ravel()
                phi_Gp = np.ascontiguousarray(phi_aGp[a], complex)
                gemv(1.0, phi_Gp, P_p, 1.0, rho_G)

            if optical_limit:
                if n==m:
                    rho_G[0] = 1.
                elif np.abs(self.e_skn[spin2][ibzkpt2, m] - self.e_skn[spin1][ibzkpt1, n]) < 1e-5:
                    rho_G[0] = 0.
                else:
                    rho_G[0] /= (self.enoshift_skn[spin2][ibzkpt2, m] - self.enoshift_skn[spin1][ibzkpt1, n])

            return rho_G

    def get_P_ai(self, k, n, spin=0, Ptmp_ai=None):

        calc = self.calc
        kd = self.calc.wfs.kd
        spos_ac = self.spos_ac
        
        ibzkpt = kd.bz2ibz_k[k]
        u = ibzkpt + kd.nibzkpts * spin
        kpt = calc.wfs.kpt_u[u]
        s = kd.sym_k[k]
        time_reversal = kd.time_reversal_k[k]
        P_ai = {}
        for a, id in enumerate(calc.wfs.setups.id_a):
            b = kd.symmetry.a_sa[s, a]
            S_c = (np.dot(spos_ac[a], kd.symmetry.op_scc[s]) - kd.symmetry.ft_sc[s] - spos_ac[b])

            #print abs(S_c.round() - S_c).max()
            #print 'S_c', abs(S_c).max()
            assert abs(S_c.round() - S_c).max() < 1e-8 ##############
            k_c = kd.ibzk_kc[kpt.k]
        
            x = np.exp(2j * pi * np.dot(k_c, S_c))
            if Ptmp_ai is None:
                P_i = np.dot(calc.wfs.setups[a].R_sii[s], kpt.P_ani[b][n]) * x
            else:
                P_i = np.dot(calc.wfs.setups[a].R_sii[s], Ptmp_ai[b]) * x
            if time_reversal:
                P_i = P_i.conj()
            P_ai[a] = P_i
        return P_ai
