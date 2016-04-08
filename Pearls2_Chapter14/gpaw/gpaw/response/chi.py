from __future__ import print_function
import sys
from time import time, ctime
import numpy as np
from math import sqrt, pi
from ase.units import Hartree, Bohr
from gpaw import extra_parameters
from gpaw.utilities.blas import gemv, scal, axpy, czher
from gpaw.mpi import world, rank, size, serial_comm
from gpaw.fd_operators import Gradient
from gpaw.response.math_func import hilbert_transform
from gpaw.response.parallel import set_communicator, \
     parallel_partition, parallel_partition_list, SliceAlongFrequency, SliceAlongOrbitals
from gpaw.response.kernel import calculate_Kxc, calculate_Kc
from gpaw.response.kernel import CoulombKernel
from gpaw.utilities.memory import maxrss
from gpaw.response.base import BASECHI

class CHI(BASECHI):
    """This class is a calculator for the linear density response function.

    Parameters:

        nband: int
            Number of bands.
        wmax: floadent
            Maximum energy for spectrum.
        dw: float
            Frequency interval.
        wlist: tuple
            Frequency points.
        q: ndarray
            Momentum transfer in reduced coordinate.
        Ecut: ndarray
            Planewave cutoff energy.
        eta: float
            Spectrum broadening factor.
        sigma: float
            Width for delta function.
    """

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
                 vcut=None,
                 ftol=1e-5,
                 txt=None,
                 xc='ALDA',
                 hilbert_trans=True,
                 time_ordered=False,
                 optical_limit=False,
                 comm=None,
                 kcommsize=None):

        BASECHI.__init__(self, calc=calc, nbands=nbands, w=w, q=q,
                         eshift=eshift, ecut=ecut,
                         density_cut=density_cut, G_plus_q=G_plus_q, eta=eta,
                         rpad=rpad, ftol=ftol, txt=txt,
                         optical_limit=optical_limit)
        
        #if vcut is None:
        #    vcut = '%dD' % self.calc.wfs.gd.pbc_c.sum()
        self.vcut = vcut

        self.xc = xc
        self.hilbert_trans = hilbert_trans
        self.full_hilbert_trans = time_ordered
        self.kcommsize = kcommsize
        self.comm = comm
        if self.comm is None:
            self.comm = world
        self.chi0_wGG = None
        
    def initialize(self, simple_version=False):

        self.printtxt('')
        self.printtxt('-----------------------------------------')
        self.printtxt('Response function calculation started at:')
        self.starttime = time()
        self.printtxt(ctime())

        BASECHI.initialize(self)

        # Frequency init
        self.dw = None
        if len(self.w_w) == 1:
            self.hilbert_trans = False

        if self.hilbert_trans:
            self.dw = self.w_w[1] - self.w_w[0]
#            assert ((self.w_w[1:] - self.w_w[:-1] - self.dw) < 1e-10).all() # make sure its linear w grid
            assert self.w_w.max() == self.w_w[-1]
            
            self.dw /= Hartree
            self.w_w /= Hartree
            self.wmax = self.w_w[-1] 
            self.wcut = self.wmax + 5. / Hartree
#            self.Nw  = int(self.wmax / self.dw) + 1
            self.Nw = len(self.w_w)
            self.NwS = int(self.wcut / self.dw) + 1
        else:
            self.Nw = len(self.w_w)
            self.NwS = 0
            if len(self.w_w) > 2:
                self.dw = self.w_w[1] - self.w_w[0]
                assert ((self.w_w[1:] - self.w_w[:-1] - self.dw) < 1e-10).all()
                self.dw /= Hartree

        self.nvalbands = self.nbands
        tmpn = np.zeros(self.nspins, dtype=int)
        for spin in range(self.nspins):
            for n in range(self.nbands):
                if (self.f_skn[spin][:, n] - self.ftol < 0).all():
                    tmpn[spin] = n
                    break
        if tmpn.max() > 0:
            self.nvalbands = tmpn.max()

        # Parallelization initialize
        self.parallel_init()

        # Printing calculation information
        self.print_chi()

        if extra_parameters.get('df_dry_run'):
            raise SystemExit

        calc = self.calc

        # For LCAO wfs
        if calc.input_parameters['mode'] == 'lcao':
            calc.initialize_positions()        
        self.printtxt('     Max mem sofar   : %f M / cpu' %(maxrss() / 1024**2))

        if simple_version is True:
            return
        # PAW part init
        # calculate <phi_i | e**(-i(q+G).r) | phi_j>
        # G != 0 part
        self.phi_aGp, self.phiG0_avp = self.get_phi_aGp(alldir=True)
        self.printtxt('Finished phi_aGp !')
        mem = np.array([self.phi_aGp[i].size * 16 /1024.**2 for i in range(len(self.phi_aGp))])
        self.printtxt('     Phi_aGp         : %f M / cpu' %(mem.sum()))

        # Calculate ALDA kernel (not used in chi0)
        R_av = calc.atoms.positions / Bohr
        if self.xc == 'RPA': #type(self.w_w[0]) is float:
            self.Kc_GG = None
            self.printtxt('RPA calculation.')
        elif self.xc == 'ALDA' or self.xc == 'ALDA_X':
            #self.Kc_GG = calculate_Kc(self.q_c,
            #                          self.Gvec_Gc,
            #                          self.acell_cv,
            #                          self.bcell_cv,
            #                          self.calc.atoms.pbc,
            #                          self.vcut)
            # Initialize a CoulombKernel instance
            kernel = CoulombKernel(vcut=self.vcut,
                                   pbc=self.calc.atoms.pbc,
                                   cell=self.acell_cv)
            self.Kc_GG = kernel.calculate_Kc(self.q_c,
                                             self.Gvec_Gc,
                                             self.bcell_cv)
            
            self.Kxc_sGG = calculate_Kxc(self.gd, # global grid
                                         self.gd.zero_pad(calc.density.nt_sG),
                                         self.npw, self.Gvec_Gc,
                                         self.gd.N_c, self.vol,
                                         self.bcell_cv, R_av,
                                         calc.wfs.setups,
                                         calc.density.D_asp,
                                         functional=self.xc,
                                         density_cut=self.density_cut)
            
            self.printtxt('Finished %s kernel ! ' % self.xc)
                
        return


    def pawstuff(self, psit_g, k, n, spin, u, ibzkpt):
        if not self.pwmode:
            if self.calc.wfs.world.size > 1 or self.kd.nbzkpts == 1:
                P_ai = self.pt.dict()
                self.pt.integrate(psit_g, P_ai, k)
            else:
                P_ai = self.get_P_ai(k, n, spin)
        else:
            # first calculate P_ai at ibzkpt, then rotate to k
            Ptmp_ai = self.pt.dict()
            kpt = self.calc.wfs.kpt_u[u]
            self.pt.integrate(kpt.psit_nG[n], Ptmp_ai, ibzkpt)
            P_ai = self.get_P_ai(k, n, spin, Ptmp_ai)
        return P_ai


    def calculate(self, seperate_spin=None):
        """Calculate the non-interacting density response function. """
        calc = self.calc
        kd = self.kd
        gd = self.gd
        sdisp_cd = gd.sdisp_cd
        ibzk_kc = kd.ibzk_kc
        bzk_kc = kd.bzk_kc
        kq_k = self.kq_k
        f_skn = self.f_skn
        e_skn = self.e_skn

        # Matrix init
        chi0_wGG = np.zeros((self.Nw_local, self.npw, self.npw), dtype=complex)
        if self.hilbert_trans:
            specfunc_wGG = np.zeros((self.NwS_local, self.npw, self.npw), dtype = complex)

        # Prepare for the derivative of pseudo-wavefunction
        if self.optical_limit:
            d_c = [Gradient(gd, i, n=4, dtype=complex).apply for i in range(3)]
            dpsit_g = gd.empty(dtype=complex)
            tmp = np.zeros((3), dtype=complex)
            rhoG0_v = np.zeros(3, dtype=complex)

            self.chi0G0_wGv = np.zeros((self.Nw_local, self.npw, 3), dtype=complex)
            self.chi00G_wGv = np.zeros((self.Nw_local, self.npw, 3), dtype=complex)

            specfuncG0_wGv = np.zeros((self.NwS_local, self.npw, 3), dtype=complex)
            specfunc0G_wGv = np.zeros((self.NwS_local, self.npw, 3), dtype=complex)
            
        use_zher = False
        if self.eta < 1e-5:
            use_zher = True

        rho_G = np.zeros(self.npw, dtype=complex)
        t0 = time()

        if seperate_spin is None:
            spinlist = np.arange(self.nspins)
        else:
            spinlist = [seperate_spin]
        
        for spin in spinlist:
            if not (f_skn[spin] > self.ftol).any():
                self.chi0_wGG = chi0_wGG
                continue

            for k in range(self.kstart, self.kend):
                k_pad = False
                if k >= self.kd.nbzkpts:
                    k = 0
                    k_pad = True
    
                # Find corresponding kpoint in IBZ
                ibzkpt1 = kd.bz2ibz_k[k]
                if self.optical_limit:
                    ibzkpt2 = ibzkpt1
                else:
                    ibzkpt2 = kd.bz2ibz_k[kq_k[k]]
    
                if self.pwmode:
                    N_c = self.gd.N_c
                    k_c = self.kd.ibzk_kc[ibzkpt1]
                    eikr1_R = np.exp(2j * pi * np.dot(np.indices(N_c).T, k_c / N_c).T)
                    k_c = self.kd.ibzk_kc[ibzkpt2]
                    eikr2_R = np.exp(2j * pi * np.dot(np.indices(N_c).T, k_c / N_c).T)
                    
                index1_g, phase1_g = kd.get_transform_wavefunction_index(self.gd.N_c - (self.pbc == False), k)
                index2_g, phase2_g = kd.get_transform_wavefunction_index(self.gd.N_c - (self.pbc == False), kq_k[k])
                
                for n in range(self.nvalbands):
                    if self.calc.wfs.world.size == 1:
                        if (self.f_skn[spin][ibzkpt1, n] - self.ftol < 0):
                            continue

                    t1 = time()
                    if self.pwmode:
                        u = self.kd.get_rank_and_index(spin, ibzkpt1)[1]
                        psitold_g = calc.wfs._get_wave_function_array(u, n, realspace=True, phase=eikr1_R)
                    else:
                        u = None
                        psitold_g = self.get_wavefunction(ibzkpt1, n, True, spin=spin)
    
                    psit1new_g = kd.transform_wave_function(psitold_g,k,index1_g,phase1_g)

                    P1_ai = self.pawstuff(psit1new_g, k, n, spin, u, ibzkpt1)

                    psit1_g = psit1new_g.conj() * self.expqr_g
    
                    for m in self.mlist:
                        if self.nbands > 1000 and m % 200 == 0:
                            print('    ', k, n, m, time() - t0, file=self.txt)
                    
                        check_focc = (f_skn[spin][ibzkpt1, n] - f_skn[spin][ibzkpt2, m]) > self.ftol
    
                        if not self.pwmode:
                            psitold_g = self.get_wavefunction(ibzkpt2, m, check_focc, spin=spin)
    
                        if check_focc:                            
                            if self.pwmode:
                                u = self.kd.get_rank_and_index(spin, ibzkpt2)[1]
                                psitold_g = calc.wfs._get_wave_function_array(u, m, realspace=True, phase=eikr2_R)
    
                            psit2_g = kd.transform_wave_function(psitold_g, kq_k[k], index2_g, phase2_g)
    
                            # zero padding is included through the FFT
                            rho_g = np.fft.fftn(psit2_g * psit1_g, s=self.nGrpad) * self.vol / self.nG0rpad
                            # Here, planewave cutoff is applied
                            rho_G = rho_g.ravel()[self.Gindex_G]
    
                            if self.optical_limit:
                                phase_cd = np.exp(2j * pi * sdisp_cd * kd.bzk_kc[kq_k[k], :, np.newaxis])
                                for ix in range(3):
                                    d_c[ix](psit2_g, dpsit_g, phase_cd)
                                    tmp[ix] = gd.integrate(psit1_g * dpsit_g)
                                rho_G[0] = -1j * np.dot(self.qq_v, tmp)

                                for ix in range(3):
                                    q2_c = np.diag((1,1,1))[ix] * self.qopt
                                    qq2_v = np.dot(q2_c, self.bcell_cv) # summation over c
                                    rhoG0_v[ix] = -1j * np.dot(qq2_v, tmp)

                            P2_ai = self.pawstuff(psit2_g, kq_k[k], m, spin, u, ibzkpt2)

                            for a, id in enumerate(calc.wfs.setups.id_a):
                                P_p = np.outer(P1_ai[a].conj(), P2_ai[a]).ravel()
                                gemv(1.0, self.phi_aGp[a], P_p, 1.0, rho_G)

                                if self.optical_limit:
                                    gemv(1.0, self.phiG0_avp[a], P_p, 1.0, rhoG0_v)

                            if self.optical_limit:
                                if np.abs(self.enoshift_skn[spin][ibzkpt2, m] -
                                          self.enoshift_skn[spin][ibzkpt1, n]) > 0.1/Hartree:
                                    rho_G[0] /= self.enoshift_skn[spin][ibzkpt2, m] \
                                                - self.enoshift_skn[spin][ibzkpt1, n]
                                    rhoG0_v /= self.enoshift_skn[spin][ibzkpt2, m] \
                                                - self.enoshift_skn[spin][ibzkpt1, n]
                                else:
                                    rho_G[0] = 0.
                                    rhoG0_v[:] = 0.
    
                            if k_pad:
                                rho_G[:] = 0.

                            if self.optical_limit:
                                rho0G_Gv = np.outer(rho_G.conj(), rhoG0_v)
                                rhoG0_Gv = np.outer(rho_G, rhoG0_v.conj())
                                rho0G_Gv[0,:] = rhoG0_v * rhoG0_v.conj()
                                rhoG0_Gv[0,:] = rhoG0_v * rhoG0_v.conj()

                            if not self.hilbert_trans:
                                if not use_zher:
                                    rho_GG = np.outer(rho_G, rho_G.conj())

                                for iw in range(self.Nw_local):
                                    w = self.w_w[iw + self.wstart] / Hartree
                                    coef = ( 1. / (w + e_skn[spin][ibzkpt1, n] - e_skn[spin][ibzkpt2, m]
                                                   + 1j * self.eta) 
                                           - 1. / (w - e_skn[spin][ibzkpt1, n] + e_skn[spin][ibzkpt2, m]
                                                   + 1j * self.eta) )
                                    C =  (f_skn[spin][ibzkpt1, n] - f_skn[spin][ibzkpt2, m]) * coef

                                    if use_zher:
                                        czher(C.real, rho_G.conj(), chi0_wGG[iw])
                                    else:
                                        axpy(C, rho_GG, chi0_wGG[iw])
                                        
                                        if self.optical_limit:
                                            axpy(C, rho0G_Gv, self.chi00G_wGv[iw])
                                            axpy(C, rhoG0_Gv, self.chi0G0_wGv[iw])
    
                            else:
                                rho_GG = np.outer(rho_G, rho_G.conj())
                                focc = f_skn[spin][ibzkpt1,n] - f_skn[spin][ibzkpt2,m]
                                w0 = e_skn[spin][ibzkpt2,m] - e_skn[spin][ibzkpt1,n]
                                scal(focc, rho_GG)
                                if self.optical_limit:
                                    scal(focc, rhoG0_Gv)
                                    scal(focc, rho0G_Gv)
    
                                # calculate delta function
                                w0_id = int(w0 / self.dw)
                                if w0_id + 1 < self.NwS:
                                    # rely on the self.NwS_local is equal in each node!
                                    if self.wScomm.rank == w0_id // self.NwS_local:
                                        alpha = (w0_id + 1 - w0/self.dw) / self.dw
                                        axpy(alpha, rho_GG, specfunc_wGG[w0_id % self.NwS_local] )

                                        if self.optical_limit:
                                            axpy(alpha, rho0G_Gv, specfunc0G_wGv[w0_id % self.NwS_local] )
                                            axpy(alpha, rhoG0_Gv, specfuncG0_wGv[w0_id % self.NwS_local] )
    
                                    if self.wScomm.rank == (w0_id+1) // self.NwS_local:
                                        alpha =  (w0 / self.dw - w0_id) / self.dw
                                        axpy(alpha, rho_GG, specfunc_wGG[(w0_id+1) % self.NwS_local] )

                                        if self.optical_limit:
                                            axpy(alpha, rho0G_Gv, specfunc0G_wGv[(w0_id+1) % self.NwS_local] )
                                            axpy(alpha, rhoG0_Gv, specfuncG0_wGv[(w0_id+1) % self.NwS_local] )

    #                            deltaw = delta_function(w0, self.dw, self.NwS, self.sigma)
    #                            for wi in range(self.NwS_local):
    #                                if deltaw[wi + self.wS1] > 1e-8:
    #                                    specfunc_wGG[wi] += tmp_GG * deltaw[wi + self.wS1]
                    if self.kd.nbzkpts == 1:
                        if n == 0:
                            dt = time() - t0
                            totaltime = dt * self.nvalbands * self.nspins
                            self.printtxt('Finished n 0 in %d seconds, estimate %d seconds left.' %(dt, totaltime) )
                        if rank == 0 and self.nvalbands // 5 > 0:
                            if n > 0 and n % (self.nvalbands // 5) == 0:
                                dt = time() - t0
                                self.printtxt('Finished n %d in %d seconds, estimate %d seconds left.'%(n, dt, totaltime-dt))
                if calc.wfs.world.size != 1:
                    self.kcomm.barrier()            
                if k == 0:
                    dt = time() - t0
                    totaltime = dt * self.nkpt_local * self.nspins
                    self.printtxt('Finished k 0 in %d seconds, estimate %d seconds left.' %(dt, totaltime))
                    
                if rank == 0 and self.nkpt_local // 5 > 0:            
                    if k > 0 and k % (self.nkpt_local // 5) == 0:
                        dt =  time() - t0
                        self.printtxt('Finished k %d in %d seconds, estimate %d seconds left.  '%(k, dt, totaltime - dt) )
        self.printtxt('Finished summation over k')

        self.kcomm.barrier()
        
        # Hilbert Transform
        if not self.hilbert_trans:
            for iw in range(self.Nw_local):
                self.kcomm.sum(chi0_wGG[iw])
                if self.optical_limit:
                    self.kcomm.sum(self.chi0G0_wGv[iw])
                    self.kcomm.sum(self.chi00G_wGv[iw])

            if use_zher:
                assert (np.abs(chi0_wGG[0,1:,0]) < 1e-10).all()
                for iw in range(self.Nw_local):
                    chi0_wGG[iw] += chi0_wGG[iw].conj().T
                    for iG in range(self.npw):
                        chi0_wGG[iw, iG, iG] /= 2.
                        assert np.abs(np.imag(chi0_wGG[iw, iG, iG])) < 1e-10 
        else:
            for iw in range(self.NwS_local):
                self.kcomm.sum(specfunc_wGG[iw])
                if self.optical_limit:
                    self.kcomm.sum(specfuncG0_wGv[iw])
                    self.kcomm.sum(specfunc0G_wGv[iw])

            if self.wScomm.size == 1:
                chi0_wGG = hilbert_transform(specfunc_wGG, self.w_w, self.Nw, self.dw, self.eta,
                                             self.full_hilbert_trans)[self.wstart:self.wend]
                self.printtxt('Finished hilbert transform !')
                del specfunc_wGG
            else:
                # redistribute specfunc_wGG to all nodes
                size = self.comm.size
                assert self.NwS % size == 0
                NwStmp1 = (rank % self.kcomm.size) * self.NwS // size
                NwStmp2 = (rank % self.kcomm.size + 1) * self.NwS // size 
                specfuncnew_wGG = specfunc_wGG[NwStmp1:NwStmp2]
                del specfunc_wGG
                
                coords = np.zeros(self.wcomm.size, dtype=int)
                nG_local = self.npw**2 // self.wcomm.size
                if self.wcomm.rank == self.wcomm.size - 1:
                    nG_local = self.npw**2 - (self.wcomm.size - 1) * nG_local
                self.wcomm.all_gather(np.array([nG_local]), coords)
        
                specfunc_Wg = SliceAlongFrequency(specfuncnew_wGG, coords, self.wcomm)
                self.printtxt('Finished Slice Along Frequency !')
                chi0_Wg = hilbert_transform(specfunc_Wg, self.w_w, self.Nw, self.dw, self.eta,
                                            self.full_hilbert_trans)[:self.Nw]
                self.printtxt('Finished hilbert transform !')
                self.comm.barrier()
                del specfunc_Wg
        
                chi0_wGG = SliceAlongOrbitals(chi0_Wg, coords, self.wcomm)
                self.printtxt('Finished Slice along orbitals !')
                self.comm.barrier()
                del chi0_Wg

                if self.optical_limit:
                    specfuncG0_WGv = np.zeros((self.NwS, self.npw, 3), dtype=complex)
                    specfunc0G_WGv = np.zeros((self.NwS, self.npw, 3), dtype=complex)
                    self.wScomm.all_gather(specfunc0G_wGv, specfunc0G_WGv)
                    self.wScomm.all_gather(specfuncG0_wGv, specfuncG0_WGv)
                    specfunc0G_wGv = specfunc0G_WGv
                    specfuncG0_wGv = specfuncG0_WGv

            if self.optical_limit:
                self.chi00G_wGv = hilbert_transform(specfunc0G_wGv, self.w_w, self.Nw, self.dw, self.eta,
                                             self.full_hilbert_trans)[self.wstart:self.wend]
                
                self.chi0G0_wGv = hilbert_transform(specfuncG0_wGv, self.w_w, self.Nw, self.dw, self.eta,
                                             self.full_hilbert_trans)[self.wstart:self.wend]

        if self.optical_limit:
            self.chi00G_wGv /= self.vol
            self.chi0G0_wGv /= self.vol

        
        self.chi0_wGG = chi0_wGG
        self.chi0_wGG /= self.vol

        self.printtxt('')
        self.printtxt('Finished chi0 !')


    def parallel_init(self):
        """Parallel initialization. By default, only use kcomm and wcomm.

        Parameters:

            kcomm:
                 kpoint communicator
            wScomm:
                 spectral function communicator
            wcomm:
                 frequency communicator
        """

        if extra_parameters.get('df_dry_run'):
            from gpaw.mpi import DryRunCommunicator
            size = extra_parameters['df_dry_run']
            world = DryRunCommunicator(size)
            rank = world.rank
            self.comm = world
        else:
            world = self.comm
            rank = self.comm.rank
            size = self.comm.size

        wcommsize = int(self.NwS * self.npw**2 * 16. / 1024**2) // 1500 # megabyte
        wcommsize += 1
        if size < wcommsize:
            raise ValueError('Number of cpus are not enough ! ')
        if self.kcommsize is None:
            self.kcommsize = world.size
        if wcommsize > size // self.kcommsize: # if matrix too large, overwrite kcommsize and distribute matrix
            self.printtxt('kcommsize is over written ! ')
            while size % wcommsize != 0:
                wcommsize += 1
            self.kcommsize = size // wcommsize
            assert self.kcommsize * wcommsize == size
            if self.kcommsize < 1:
                raise ValueError('Number of cpus are not enough ! ')

        self.kcomm, self.wScomm, self.wcomm = set_communicator(world, rank, size, self.kcommsize)

        if self.kd.nbzkpts >= world.size:
            self.nkpt_reshape = self.kd.nbzkpts
            self.nkpt_reshape, self.nkpt_local, self.kstart, self.kend = parallel_partition(
                               self.nkpt_reshape, self.kcomm.rank, self.kcomm.size, reshape=True, positive=True)
            self.mband_local = self.nvalbands
            self.mlist = np.arange(self.nbands)
        else:
            # if number of kpoints == 1, use band parallelization
            self.nkpt_local = self.kd.nbzkpts
            self.kstart = 0
            self.kend = self.kd.nbzkpts
            self.nkpt_reshape = self.kd.nbzkpts

            self.nbands, self.mband_local, self.mlist = parallel_partition_list(
                               self.nbands, self.kcomm.rank, self.kcomm.size)

        if self.NwS % size != 0:
            self.NwS -= self.NwS % size
            
        self.NwS, self.NwS_local, self.wS1, self.wS2 = parallel_partition(
                               self.NwS, self.wScomm.rank, self.wScomm.size, reshape=False)

        if self.hilbert_trans:
            self.Nw, self.Nw_local, self.wstart, self.wend =  parallel_partition(
                               self.Nw, self.wcomm.rank, self.wcomm.size, reshape=True)
        else:
            if self.Nw > 1:
#                assert self.Nw % (self.comm.size / self.kcomm.size) == 0
                self.wcomm = self.wScomm
                self.Nw, self.Nw_local, self.wstart, self.wend =  parallel_partition(
                               self.Nw, self.wcomm.rank, self.wcomm.size, reshape=False)
            else:
                # if frequency point is too few, then dont parallelize
                self.wcomm = serial_comm
                self.wstart = 0
                self.wend = self.Nw
                self.Nw_local = self.Nw

        return

    def print_chi(self):

        printtxt = self.printtxt
        printtxt('Use Hilbert Transform: %s' %(self.hilbert_trans) )
        printtxt('Calculate time-ordered Response Function: %s' %(self.full_hilbert_trans) )
        printtxt('')
        printtxt('Number of frequency points   : %d' %(self.Nw) )
        if self.hilbert_trans:
            printtxt('Number of specfunc points    : %d' % (self.NwS))
        printtxt('')
        printtxt('Parallelization scheme:')
        printtxt('     Total cpus      : %d' %(self.comm.size))
        if self.kd.nbzkpts == 1:
            printtxt('     nbands parsize  : %d' %(self.kcomm.size))
        else:
            printtxt('     kpoint parsize  : %d' %(self.kcomm.size))
            if self.nkpt_reshape > self.kd.nbzkpts:
                self.printtxt('        kpoints (%d-%d) are padded with zeros' % (self.kd.nbzkpts, self.nkpt_reshape))

        if self.hilbert_trans:
            printtxt('     specfunc parsize: %d' %(self.wScomm.size))
        printtxt('     w parsize       : %d' %(self.wcomm.size))
        printtxt('')
        printtxt('Memory usage estimation:')
        printtxt('     chi0_wGG        : %f M / cpu' %(self.Nw_local * self.npw**2 * 16. / 1024**2) )
        if self.hilbert_trans:
            printtxt('     specfunc_wGG    : %f M / cpu' %(self.NwS_local *self.npw**2 * 16. / 1024**2) )

